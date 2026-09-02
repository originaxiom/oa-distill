import json, os, sys
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
def say(*a): print(*a); sys.stdout.flush()
def write(name, data):
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, name + '.json'), 'w') as f: json.dump(data, f, indent=1, sort_keys=True, default=str)
    say(f"[{name}] written")
