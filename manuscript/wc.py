import re, sys
def wc(t):
    t = re.sub(r"\$\$.*?\$\$", " ", t, flags=re.S); t = re.sub(r"\$[^$]*\$", " X ", t)
    t = re.sub(r"\{\{[^}]*\}\}", " ", t); t = re.sub(r"\|.*\|", " ", t); t = re.sub(r":::[^\n]*", " ", t)
    return len(re.findall(r"[A-Za-z0-9][\w\-–’']*", t))
tot = 0
for f in sys.argv[1:]:
    s = open(f).read()
    if "part4" in f: s = s.split("**Author Contributions")[0]
    if "part1" in f: s = s[s.index("**Abstract"):]
    n = wc(s); tot += n; print(f"{n:6d} {f}")
    if "-v" in sys.argv[0]: pass
print("TOTAL", tot)
