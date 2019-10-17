#!/usr/bin/env python
"""
Script to rename pseudodojo files under revision control (git)
"""
import sys
import os
import json


def main():
    usage = "Usage: dojomv.py old_pseudo new_pseudo"
    if "--help" in sys.argv or "-h" in sys.argv:
        print(usage)
        return 1
    try:
        old, new = sys.argv[1:]
    except:
        print(usage)
        return 1

    old = os.path.abspath(old)
    old_root, old_ext = os.path.splitext(old)

    new = os.path.abspath(new)
    new_root, new_ext = os.path.splitext(new)

    assert old_ext == new_ext and old_ext == ".psp8"
    exts = ["djrepo", "out", "psp8", "in"]
    # Git-remove new files (if present)
    for ext in exts:
        path = os.path.join(new_root, "." + ext)
        if os.path.exists(path):
            cmd = "git rm %s" % path
            ret = os.system(cmd)
            print(cmd, "[%s]" % ret)

    # Git rename old files.
    for ext in exts:
        old_path = old_root + "." + ext
        new_path = new_root + "." + ext
        if os.path.exists(old_path):
            cmd = "git mv %s %s" % (old_path, new_path)
            ret = os.system(cmd)
            print(cmd, "[%s]" % ret)
        else:
            print("Warning: old_path %s does not exist" % old_path)

    # Change basename in djrepo file.
    djrepo = new_root + ".djrepo"
    with open(djrepo, "r") as fh:
        d = json.load(fh)
        d["basename"] = os.path.basename(new_root + ".psp8")

    with open(djrepo, "w") as fh:
        json.dump(d, fh, indent=-1, sort_keys=True) #, cls=MontyEncoder)

    return 0

if __name__ == "__main__":
    sys.exit(main())
