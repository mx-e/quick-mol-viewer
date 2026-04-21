# quick-mol-viewer

A single-page, static molecule viewer. All state lives in the URL — no backend, no uploads, no build step. Paste a link, render a molecule.

**Live:** https://mx-e.github.io/quick-mol-viewer/

## What it does

Three content modes, all served from one HTML file:

- **Single snapshot** — one molecule, one frame.
- **Trajectory** — one molecule, many frames. Play/pause, scrub, variable speed.
- **Collection** — multiple independent molecules (different atom counts OK). Navigate one-at-a-time or view a 3×3 grid (press `G`).

## Usage

Two ways to load a molecule:

1. **URL fragment** — `https://mx-e.github.io/quick-mol-viewer/#d=<base64url-payload>`
2. **Paste field** — click `paste` (top-right) or press `P`, paste either the full URL or just the payload.

The payload is a base64url-encoded, optionally gzip-compressed, little-endian binary blob: header (version, flags, atom/frame counts, scale) followed by atomic numbers and quantised int16 coordinates. Delta-encoded trajectory frames and heterogeneous collections are supported. Copy-paste artifacts (line wraps, backslashes, `%XX`, smart quotes) are stripped automatically before decoding.

Full byte layout and reference decoder: see the [build spec](./SPEC.md) if included, or the `decodeMolPayload` function in `index.html`.

## Keyboard

| key | action |
| --- | --- |
| `Space` | play/pause trajectory |
| `←` / `→` | step frame (trajectory) / molecule (collection, single) / page of 9 (collection, grid) |
| `G` | toggle grid view (collection) |
| `P` | open paste dialog |
| `Esc` | close paste dialog |
| `⌘/Ctrl+Enter` | submit paste |

## Encoder

The Python encoder (`mol_url.py`) lives separately. A matching JS encoder is trivial to lift from the test files in `/tmp/mol_*_test.mjs` if you want to generate URLs client-side.

## Credits

Rendering by [3Dmol.js](https://3dmol.csb.pitt.edu/).
