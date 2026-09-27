r"""Renders the trailer's new scenes (see ``scenes.py``) as stills and as animations.

For each scene: ``build/stills/<name>.png`` (the finished frame), ``build/frames/<name>/f*.png``
(1280x720, 30 fps, over the shot's length at 96 bpm; the last 12 % holds the finished frame), and
a clip in ``build/clips/``.  Run ``data.py`` first.  One scene at a time on one core::

    nice -n 19 python tools/trailer/render.py [--stills-only] [scene ...]

The clip is H.264 in MP4 from the PNG frames when the ffmpeg on the PATH (or in $FFMPEG) can
do that.  Otherwise, as with the ffmpeg bundled with Playwright, which reads only piped JPEG and
writes only VP8, the frames go through quality-95 JPEG into a WebM for review; the PNG frames
remain the masters either way.
"""
import os
import shutil
import subprocess
import sys

from common import BUILD, BG, FPS, SPB, setup_matplotlib

plt = setup_matplotlib()
import scenes  # noqa: E402  (needs the fonts registered first)

HOLD = 0.12
PLAYWRIGHT_FFMPEG = '/opt/pw-browsers/ffmpeg-1011/ffmpeg-linux'


def ffmpeg():
    for exe in (os.environ.get('FFMPEG'), shutil.which('ffmpeg'), PLAYWRIGHT_FFMPEG):
        if exe and os.path.exists(exe):
            return exe
    raise SystemExit('no ffmpeg found: set $FFMPEG')


def frame(draw, u, path):
    fig = plt.figure(figsize=(16, 9), dpi=80, facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis('off')
    draw(fig, ax, u)
    fig.savefig(path, facecolor=BG)
    plt.close(fig)


def encode(frames, n, name):
    exe = ffmpeg()
    caps = subprocess.run([exe, '-hide_banner', '-encoders'], capture_output=True, text=True).stdout
    demux = subprocess.run([exe, '-hide_banner', '-demuxers'], capture_output=True, text=True).stdout
    clips = BUILD / 'clips'
    clips.mkdir(exist_ok=True)
    if 'libx264' in caps and ' image2 ' in demux:
        dest = clips / (name + '.mp4')
        subprocess.run([exe, '-y', '-loglevel', 'error', '-framerate', str(FPS), '-i', str(frames / 'f%04d.png'),
                        '-c:v', 'libx264', '-crf', '16', '-pix_fmt', 'yuv420p', str(dest)], check=True)
        return dest
    from PIL import Image
    cat = frames.with_suffix('.mjpeg')
    with open(cat, 'wb') as fh:
        for i in range(n):
            Image.open(frames / ('f%04d.png' % i)).convert('RGB').save(fh, 'JPEG', quality=95)
    dest = clips / (name + '.webm')
    subprocess.run([exe, '-y', '-loglevel', 'error', '-f', 'image2pipe', '-c:v', 'mjpeg', '-framerate', str(FPS),
                    '-i', 'file:' + str(cat), '-c:v', 'libvpx', '-b:v', '3M', '-auto-alt-ref', '0',
                    '-pix_fmt', 'yuv420p', str(dest)], check=True)
    cat.unlink()
    return dest


def render(name, stills_only=False):
    draw, beats, _ = scenes.SCENES[name]
    (BUILD / 'stills').mkdir(parents=True, exist_ok=True)
    frame(draw, 1.0, BUILD / 'stills' / (name + '.png'))
    if stills_only:
        return
    frames = BUILD / 'frames' / name
    frames.mkdir(parents=True, exist_ok=True)
    n = int(round(beats * SPB * FPS))
    for i in range(n):
        frame(draw, min(1.0, i / ((n - 1) * (1 - HOLD))), frames / ('f%04d.png' % i))
    dest = encode(frames, n, name)
    print('%s: %d frames, %s (%d KB)' % (name, n, dest.name, dest.stat().st_size // 1024), flush=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    stills_only = '--stills-only' in args
    for name in [a for a in args if not a.startswith('--')] or list(scenes.SCENES):
        render(name, stills_only)
