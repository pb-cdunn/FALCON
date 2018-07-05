from __future__ import absolute_import
from __future__ import division

from pypeflow.io import (
        syscall, capture, cd,
        mkdirs, symlink, rm, touch, filesize, exists_and_not_empty) # needed?
import contextlib
import io
import logging
import os
import pprint
import sys

if sys.version_info >= (3, 0):
    NativeIO = io.StringIO
else:
    NativeIO = io.BytesIO

LOG = logging.getLogger()


def log(*msgs):
    LOG.debug(' '.join(repr(m) for m in msgs))


def eng(number):
    return '{:.1f}MB'.format(number / 2**20)


class Percenter(object):
    """Report progress by golden exponential.

    Usage:
        counter = Percenter('mystruct', total_len(mystruct))

        for rec in mystruct:
            counter(len(rec))
    """
    def __init__(self, name, total, log=LOG.info, units='units'):
        if sys.maxint == total:
            log('Counting {} from "{}"'.format(units, name))
        else:
            log('Counting {:,d} {} from\n  "{}"'.format(total, units, name))
        self.total = total
        self.log = log
        self.name = name
        self.units = units
        self.call = 0
        self.count = 0
        self.next_count = 0
        self.a = 1 # double each time
    def __call__(self, more, label=''):
        self.call += 1
        self.count += more
        if self.next_count <= self.count:
            self.a = 2 * self.a
            self.a = max(self.a, more)
            self.a = min(self.a, (self.total-self.count), round(self.total/10.0))
            self.next_count = self.count + self.a
            if self.total == sys.maxint:
                msg = '{:>10} count={:15,d} {}'.format(
                    '#{:,d}'.format(self.call), self.count, label)
            else:
                msg = '{:>10} count={:15,d} {:6.02f}% {}'.format(
                    '#{:,d}'.format(self.call), self.count, 100.0*self.count/self.total, label)
            self.log(msg)
    def finish(self):
        self.log('Counted {:,d} {} in {} calls from:\n  "{}"'.format(
            self.count, self.units, self.call, self.name))


def FilePercenter(fn, log=LOG.info):
    if '-' == fn or not fn:
        size = sys.maxint
    else:
        size = filesize(fn)
    return Percenter(fn, size, log, units='bytes')

@contextlib.contextmanager
def open_progress(fn, mode='r', log=LOG.info):
    """
    Usage:
        with open_progress('foo', log=LOG.info) as stream:
            for line in stream:
                use(line)

    That will log progress lines.
    """
    def get_iter(stream, progress):
        for line in stream:
            progress(len(line))
            yield line

    fp = FilePercenter(fn, log=log)
    with open(fn, mode=mode) as stream:
        yield get_iter(stream, fp)
    fp.finish()

def read_as_msgpack(stream):
    import msgpack
    content = stream.read()
    log('  Read {} as msgpack'.format(eng(len(content))))
    return msgpack.loads(content)


def read_as_json(stream):
    import json
    content = stream.read()
    log('  Read {} as json'.format(eng(len(content))))
    return json.loads(content)


def write_as_msgpack(stream, val):
    import msgpack
    content = msgpack.dumps(val)
    log('  Serialized to {} as msgpack'.format(eng(len(content))))
    stream.write(content)


def write_as_json(stream, val):
    import json
    content = json.dumps(val, indent=2, separators=(',', ': '))
    log('  Serialized to {} as json'.format(eng(len(content))))
    stream.write(content)


def deserialize(fn):
    log('Deserializing from {!r}'.format(fn))
    with open(fn) as ifs:
        log('  Opened for read: {!r}'.format(fn))
        if fn.endswith('.msgpack'):
            val = read_as_msgpack(ifs)
        elif fn.endswith('.json'):
            val = read_as_json(ifs)
        else:
            raise Exception('Unknown extension for {!r}'.format(fn))
    log('  Deserialized {} records'.format(len(val)))
    return val


def serialize(fn, val):
    """Assume dirname exists.
    """
    log('Serializing {} records'.format(len(val)))
    mkdirs(os.path.dirname(fn))
    with open(fn, 'w') as ofs:
        log('  Opened for write: {!r}'.format(fn))
        if fn.endswith('.msgpack'):
            write_as_msgpack(ofs, val)
        elif fn.endswith('.json'):
            write_as_json(ofs, val)
            ofs.write('\n') # for vim
        else:
            raise Exception('Unknown extension for {!r}'.format(fn))


def yield_abspath_from_fofn(fofn_fn):
    """Yield each filename.
    Relative paths are resolved from the FOFN directory.
    'fofn_fn' can be .fofn, .json, .msgpack
    """
    try:
        fns = deserialize(fofn_fn)
    except:
        #LOG('las fofn {!r} does not seem to be JSON; try to switch, so we can detect truncated files.'.format(fofn_fn))
        fns = open(fofn_fn).read().strip().split()
    try:
        basedir = os.path.dirname(fofn_fn)
        for fn in fns:
            if not os.path.isabs(fn):
                fn = os.path.abspath(os.path.join(basedir, fn))
            yield fn
    except Exception:
        LOG.error('Problem resolving paths in FOFN {!r}'.format(fofn_fn))
        raise


def rmdirs(*dirnames):
    for d in dirnames:
        assert os.path.normpath(d.strip()) not in ['.', '', '/']
    syscall('rm -rf {}'.format(' '.join(dirnames)))

def rmdir(d):
    rmdirs(d)

def rm_force(*fns):
    for fn in fns:
        if os.path.exists(fn):
            os.unlink(fn)
