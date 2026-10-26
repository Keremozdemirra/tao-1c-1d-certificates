"""Convert a certificate JSON to the plain-text input read by the *_verify_direct.c programs.

Usage: python3 -I to_txt.py certificate.json out.txt
"""
import json
import sys

d = json.load(open(sys.argv[1]))
if "v_int" in d:
    head, vals = "%d" % len(d["v_int"]), d["v_int"]
elif "weights_int" in d:
    head, vals = "%d %d" % (d["K"], len(d["weights_int"])), d["weights_int"]
else:
    head, vals = "%d %d" % (len(d["a_int"]), d["D"]), d["a_int"]
open(sys.argv[2], "w").write(head + "\n" + "\n".join(map(str, vals)) + "\n")
