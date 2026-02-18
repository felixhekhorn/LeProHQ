"""Parse `numbers-for-felix.txt` to raw data."""

import re
import pathlib

# read
here = pathlib.Path(__file__).parent
p = here / "numbers-for-felix.txt"
cnt = p.read_text().splitlines()
f2 = cnt[3:34:3]
fl = cnt[37::3]

# parse
f2p = re.compile(r"([+-]\.\d+(?:e[+-]\d)?)([+-]\.\d+(?:e[+-]\d)?)")
raw_f2 = []
for l in f2:
    raw_f2.append(f2p.match(l).groups())

flp = re.compile(r"([+-]\.\d+(?:e[+-]\d)?).+([+-]\.\d+(?:e[+-]\d)?)")
raw_fl = []
for l in fl:
    raw_fl.append(reversed(l.split("*ln(mu^2/m^2)")))

# glue back
cnt = []
for t, l in zip(raw_f2, raw_fl):
    cnt.append(" ".join([*t, *l]))

(here / "tab2-mod.txt").write_text("\n".join(cnt))
