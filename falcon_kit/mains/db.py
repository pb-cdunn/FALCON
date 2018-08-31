"""
This is meant to be used for LA4Falcon_pre/post hooks.
dbdir is probably /dev/shm.
"""
import os, shutil, sys

# We will ignore track files (.anno/.data).
suffixes = ('.idx', '.bps')

def log(msg):
    print(msg)
def rm(bn, dn):
    """Remove bn from directory.
    Skip silently if not found.
    Leave the directory tree.
    """
    fn = os.path.join(dn, bn)
    if os.path.exists(fn):
        log('rm -f "{}"'.format(fn))
        os.remove(fn)
def cp(bn, src_dn, dst_dn):
    """Copy bn from src to dst.
    Create dirs for dst_dn as needed.
    Over-write if exists in dst.
    Raise Exception if bn is not found in src_dn.
    """
    src_fn = os.path.join(src_dn, bn)
    dst_fn = os.path.join(dst_dn, bn)
    if not os.path.exists(src_fn):
        msg = 'Nothing found at "{}"'.format(src_fn)
        raise Exception(msg)
    if not os.path.isdir(dst_dn):
        log('mkdir -p "{}"'.format(dst_dn))
        os.makedirs(dst_dn)
    if os.path.exists(dst_fn):
        log('WARNING: {!r} already exists. Deleting and re-copying.'.format(dst_fn))
        rm(bn, dst_dn)
    log('cp -f "{}" "{}"'.format(src_fn, dst_fn))
    shutil.copy2(src_fn, dst_fn)
def clean(db, dbdir):
    """
    Remove db and dot-db files from dbdir.
    Assume the same basename was used.
    """
    bn = os.path.basename(db)
    assert bn.endswith('.db'), '{} does not end in .db'.format(bn)
    dbname = bn[:-3] # drop .db
    rm(bn, dbdir)
    for suffix in suffixes:
        bn = '.'+dbname+suffix
        rm(bn, dbdir)
def copy(db, dbdir):
    """
    Copy db and dot-db files into dbdir.
    (dbdir is probably /dev/shm.)
    """
    dn, bn = os.path.split(db)
    assert bn.endswith('.db'), '{} does not end in .db'.format(bn)
    dbname = bn[:-3] # drop .db
    cp(bn, dn, dbdir)
    for suffix in suffixes:
        bn = '.'+dbname+suffix
        cp(bn, dn, dbdir)
def main(prog, subcmd, db, dbdir):
    cmd2func = {'clean': clean, 'copy': copy}
    func = cmd2func[subcmd]
    func(db, dbdir)

if __name__ == "__main__":
    main(*sys.argv)  # pylint: disable=no-value-for-parameter
