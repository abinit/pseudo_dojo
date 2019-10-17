#!/usr/bin/env python
"""
Generate new directory for new XC by copying all input files
and by changing the XC flag in the input files.
"""
import sys
import os
import shutil


def my_copytree(src, dst, symlinks=False, ignore=None):
    """Hacked version. Based on https://docs.python.org/2/library/shutil.html"""
    names = os.listdir(src)
    if ignore is not None:
        ignored_names = ignore(src, names)
    else:
        ignored_names = set()

    os.makedirs(dst)
    errors = []
    for name in names:
        if name in ignored_names:
            continue

        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)

        # Here we select only the input files.
        if os.path.isfile(srcname) and not name.endswith(".in"):
            continue

        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                my_copytree(srcname, dstname, symlinks, ignore)
            else:
                shutil.copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error) as why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except Error as err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    #except WindowsError:
        # can't copy file access times on Windows
        #pass
    except OSError as why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise Error(errors)


def main():
    usage = "Usage: newtable.py srcdir newdir newxc"
    if "--help" in sys.argv or "-h" in sys.argv:
        print(usage)
        return 1
    try:
        src, dest, xc = sys.argv[1:]
    except:
        print(usage)
        return 1

    if os.path.exists(dest):
        raise RuntimeError("Destination %s already exists" % dest)

    print("Will copy input files from %s to %s" % (src, dest))
    my_copytree(src, dest, symlinks=False, ignore=None)

    print("Setting iexc to" , xc)
    for dirpath, dirnames, filenames in os.walk(dest):
        for f in filenames:
            if not f.endswith(".in"): continue
            path = os.path.join(dirpath, f)
            with open(path ,"r") as fh:
                lines = fh.readlines()
                """
                # atsym z nc nv iexc psfile
                Ag 47 6 4 4 psp8
                """
                assert lines[0].strip() == "# atsym z nc nv iexc psfile"
                tokens = lines[1].split()
                tokens[4] = xc
                lines[1] = " ".join(tokens) + "\n"

            with open(path, "w") as fh:
                fh.writelines(lines)

    return 0


if __name__ == "__main__":
    sys.exit(main())
