# Notebook route audit

Runs every notebook in a copy of the repository with an IPython startup hook that logs each call
where `strategy='auto'` hands a request to the ladder, so that a change to the dispatch can be
checked against all of them.  Used for #70 (record in `../issue070_auto_ladder_route/`), and with
hooks of their own for #122 and #120 (patches in `../issue122_clamp_refusal/` and
`../issue120_auto_tight/`).

- `audit_w.py MIRROR LOG PREFIXES` runs the notebooks whose names start with one of the
  comma-separated prefixes, from `MIRROR/notebooks`, with `PYTHONPATH=MIRROR/src`,
  `MAGNUS_PAPER_CACHE_ONLY=1` and `IPYTHONDIR=MIRROR/../ipy`.  It writes each executed notebook
  to `MIRROR/../audited/` and prints one status line per notebook.  Five of them with disjoint
  prefixes cover the 29 notebooks in about 15 minutes.
- `audit.py MIRROR LOG [DONE]` is the single-process version; `DONE` lists prefixes to skip.
- `ipy/profile_default/startup/00-route.py` is the #70 hook: copy it under `MIRROR/../ipy/`.
  It appends one JSON line per hand-over to `$ROUTE_LOG`, and a first line naming the module it
  hooked.

Worth rerunning for any change to `_auto_prefers_ladder` or to the engines behind it.  Check the
first line of the log before trusting an empty one: a relative path once made an audit log
nothing.
