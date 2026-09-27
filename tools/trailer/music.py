r"""The trailer's music, synthesized here from sine waves and noise: an original piece, nothing
sampled, so there is nothing to license.

96 bpm in D minor, following the music cue of each act in ``trailer.json``'s storyboard:

1. Two journeys: a low drone, and a soft tone per track (one for vacuum, one for matter).
2. The burden: a pulse enters, one hit per beat, a tick on the off-beats.
3. The reveal: the first full chord with a low boom, then the groove.
4. Just the Hamiltonian: a rising arpeggio over the groove.
5. Capabilities: the full groove, building.
6. Beyond the textbook: the climax; the melody **plays a computed probability**: the Earth
   spectrum of the "Fast" diagram (P(nu_mu -> nu_e), 1-30 GeV, from ``build/diagrams.json``),
   sampled on eighth notes and mapped onto the D-minor pentatonic scale.
7. Accurate, fast, flexible: three chord hits, resolving.
8. Get it: the final chord, D major, with a long fade.

    python tools/trailer/music.py          # -> build/music.wav, as long as the cut
"""
import json
import wave

import numpy as np
from scipy.signal import fftconvolve

from common import BUILD, HERE, SPB

SR = 44100
T = json.loads((HERE / 'trailer.json').read_text())
rng = np.random.default_rng(99)


def hz(note):
    """MIDI note number -> frequency."""
    return 440.0 * 2 ** ((note - 69) / 12)


D2, A2, D3, F3, A3, C4, D4, F4, A4, D5, A5 = 38, 45, 50, 53, 57, 60, 62, 65, 69, 74, 81
CHORDS = {'Dm': [50, 53, 57, 62], 'Bb': [46, 50, 53, 58], 'F': [45, 48, 53, 57], 'C': [48, 52, 55, 60],
          'D': [50, 54, 57, 62]}
PROG = ['Dm', 'Bb', 'F', 'C']                       # one chord per bar (4 beats)
PENTA = [0, 3, 5, 7, 10]                            # D minor pentatonic, as offsets from D


def env(n, a, r, sustain=1.0):
    """Attack-release envelope over n samples (a, r in seconds)."""
    t = np.arange(n) / SR
    e = np.minimum(1.0, t / max(a, 1e-4)) * sustain
    tail = np.clip((n / SR - t) / max(r, 1e-4), 0, 1)
    return e * tail


def tone(f, dur, a=0.01, r=0.3, harm=(1.0, 0.3, 0.1), detune=0.0):
    n = int(dur * SR)
    t = np.arange(n) / SR
    y = sum(h * np.sin(2 * np.pi * f * (k + 1) * t * (1 + detune * (k % 2))) for k, h in enumerate(harm))
    return y * env(n, a, r)


def pluck(f, dur=0.45):
    n = int(dur * SR)
    t = np.arange(n) / SR
    return (np.sin(2 * np.pi * f * t) + 0.35 * np.sin(4 * np.pi * f * t)) * np.exp(-t * 7.0) * env(n, 0.004, 0.05)


def pad(notes, dur, a=0.6, r=1.2):
    return sum(tone(hz(m), dur, a, r, harm=(1.0, 0.18, 0.05), detune=0.002) for m in notes) / len(notes)


def kick(dur=0.35):
    n = int(dur * SR)
    t = np.arange(n) / SR
    f = 45 + 90 * np.exp(-t * 30)
    return np.sin(2 * np.pi * np.cumsum(f) / SR) * np.exp(-t * 9)


def noise_hit(dur, decay, hp=True):
    n = int(dur * SR)
    y = rng.standard_normal(n)
    if hp:
        y = np.diff(y, prepend=0.0)                 # a crude high-pass: brighter, hat-like
    return y * np.exp(-np.arange(n) / SR * decay)


def boom(dur=2.5):
    n = int(dur * SR)
    t = np.arange(n) / SR
    return (np.sin(2 * np.pi * 36 * t) * np.exp(-t * 1.6) + 0.15 * noise_hit(dur, 2.5, hp=False))


class Track:
    def __init__(self, seconds):
        self.y = np.zeros((int(seconds * SR) + SR * 3, 2))

    def add(self, t, sig, gain=1.0, pan=0.0):
        i = int(t * SR)
        if i >= len(self.y):
            return
        sig = sig[:len(self.y) - i]
        self.y[i:i + len(sig), 0] += sig * gain * np.sqrt(0.5 * (1 - pan))
        self.y[i:i + len(sig), 1] += sig * gain * np.sqrt(0.5 * (1 + pan))


def acts():
    """[(act, start s, end s)] from the shots."""
    out, t = [], 0.0
    for s in T['shots']:
        d = s['beats'] * SPB
        if out and out[-1][0] == s['act']:
            out[-1][2] = t + d
        else:
            out.append([s['act'], t, t + d])
        t += d
    return out


def melody_from_probability(n_notes):
    """The Earth spectrum of the 'Fast' diagram, resampled to n_notes and mapped to pitches of the
    D-minor pentatonic scale over two octaves (low probability -> low note)."""
    d = json.loads((BUILD / 'diagrams.json').read_text())['fast']
    P = np.asarray(d['P'])
    p = np.interp(np.linspace(0, len(P) - 1, n_notes), np.arange(len(P)), P)
    steps = np.round(p / p.max() * 9).astype(int)             # 10 scale steps = 2 octaves
    return [D4 + 12 * (s // 5) + PENTA[s % 5] for s in steps]


def compose(seconds):
    tr, wet = Track(seconds), Track(seconds)                   # wet: sent to the reverb
    B = SPB
    for act, t0, t1 in acts():
        nb = int(round((t1 - t0) / B))
        beat = lambda k: t0 + k * B                            # noqa: E731
        if act == 1:
            tr.add(t0, tone(hz(D2), t1 - t0 + 1.0, a=2.5, r=2.0, harm=(1, 0.5, 0.2)), 0.10)
            tr.add(t0, tone(hz(A2), t1 - t0 + 1.0, a=4.0, r=2.0, harm=(1, 0.3)), 0.05)
            for k in range(0, nb, 2):                          # one tone per track, alternating
                wet.add(beat(k), tone(hz(A5 if k % 4 == 0 else D5), 1.2, a=0.005, r=1.1, harm=(1, 0.1)),
                        0.12, pan=-0.4 if k % 4 == 0 else 0.4)
        elif act == 2:
            tr.add(t0, tone(hz(D2), t1 - t0, a=0.5, r=0.5, harm=(1, 0.4, 0.2)), 0.10)
            for k in range(nb):
                tr.add(beat(k), kick(), 0.35)
                tr.add(beat(k + 0.5), noise_hit(0.08, 60), 0.05, pan=0.3)
        elif act == 3:
            tr.add(t0, boom(), 0.9)
            wet.add(t0, pad(CHORDS['Dm'] + [D5], (t1 - t0) + 0.5, a=0.05, r=1.5), 0.8)
            for k in range(4, nb):
                tr.add(beat(k), kick(), 0.7)
                if k % 2:
                    wet.add(beat(k), noise_hit(0.18, 18), 0.14)
                tr.add(beat(k + 0.5), noise_hit(0.06, 70), 0.05, pan=-0.3)
        elif act in (4, 5, 6):
            for k in range(nb):
                chord = CHORDS[PROG[(k // 4) % 4]]
                if k % 4 == 0:
                    wet.add(beat(k), pad(chord, 4 * B + 0.4, a=0.2, r=0.8), 0.28 + 0.06 * (act - 4))
                    tr.add(beat(k), tone(hz(chord[0] - 12), 4 * B, a=0.02, r=0.4, harm=(1, 0.5, 0.25)), 0.3)
                tr.add(beat(k), kick(), 0.6 + 0.05 * (act - 4))
                if k % 2:
                    wet.add(beat(k), noise_hit(0.18, 18), 0.14)
                for h in (0.5,) if act == 4 else (0.25, 0.5, 0.75):
                    tr.add(beat(k + h), noise_hit(0.05, 80), 0.04 + 0.01 * (act - 4), pan=0.25)
                if act == 4:                                   # rising arpeggio, sixteenths
                    for q in range(4):
                        note = chord[q % 4] + 12 * (1 + (k % 4 >= 2))
                        wet.add(beat(k + q / 4), pluck(hz(note)), 0.10 + 0.04 * k / nb, pan=0.2 * (q - 1.5))
                if act == 5 and k >= nb // 2:                  # building: eighth-note plucks join
                    for q in range(2):
                        wet.add(beat(k + q / 2), pluck(hz(chord[(k + q) % 4] + 12)), 0.08)
            if act == 6:                                       # the melody plays P(nu_mu -> nu_e)
                notes = melody_from_probability(2 * nb)
                for q, m in enumerate(notes):
                    wet.add(t0 + q * B / 2, tone(hz(m), B * 0.55, a=0.01, r=0.25, harm=(1, 0.45, 0.2, 0.1)), 0.2)
        elif act == 7:
            for k, name in zip((0, 4, 8), ('Bb', 'C', 'Dm')):
                tr.add(beat(k), boom(1.4), 0.45)
                tr.add(beat(k), kick(), 0.6)
                wet.add(beat(k), pad(CHORDS[name] + [CHORDS[name][0] + 12], 4 * B, a=0.01, r=1.6), 0.45)
        elif act == 8:
            tr.add(t0, boom(3.0), 0.6)
            wet.add(t0, pad(CHORDS['D'] + [D5, 57 + 12], t1 - t0 + 2.0, a=0.3, r=(t1 - t0) * 0.8), 0.55)
            tr.add(t0, tone(hz(D2), t1 - t0 + 2.0, a=0.3, r=(t1 - t0) * 0.8, harm=(1, 0.4)), 0.25)
    return tr.y, wet.y


def write(seconds=None):
    """Composes, mixes (with a simple reverb), masters to -1 dBFS and writes build/music.wav."""
    if seconds is None:
        seconds = sum(s['beats'] for s in T['shots']) * SPB
    dry, wet = compose(seconds)
    n = int(2.2 * SR)                                          # reverb: decaying noise, stereo
    ir = rng.standard_normal((n, 2)) * np.exp(-np.arange(n) / SR * 3.2)[:, None]
    ir /= np.sqrt(np.sum(ir**2, axis=0))
    rev = np.stack([fftconvolve(wet[:, c], ir[:, c])[:len(wet)] for c in range(2)], 1)
    mix = dry + 0.65 * wet + 0.5 * rev
    mix = mix[:int(seconds * SR)]
    fade = int(1.5 * SR)                                       # gentle ends, matching the picture
    mix[:int(0.05 * SR)] *= np.linspace(0, 1, int(0.05 * SR))[:, None]
    mix[-fade:] *= np.linspace(1, 0, fade)[:, None] ** 2
    mix = np.tanh(1.6 * mix / np.max(np.abs(mix))) / np.tanh(1.6)   # soft limiter
    mix *= 10 ** (-1 / 20)
    out = BUILD / 'music.wav'
    with wave.open(str(out), 'wb') as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes((mix * 32767).astype('<i2').tobytes())
    return out


if __name__ == '__main__':
    p = write()
    print('%s: %.1f s' % (p.relative_to(HERE.parents[1]), wave.open(str(p)).getnframes() / SR))
