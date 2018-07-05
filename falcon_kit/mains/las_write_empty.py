import argparse, logging, struct, sys

LOG = logging.getLogger()


def write_empty_las(stream):
    """Empty las means novl=0, tspace=0.
    stream should be opened for binary writing.
    """
    data = struct.pack('qi', 0, 1000)
    stream.write(data)


def run(las_fn):
    LOG.info('Writing empty las file {!r}.'.format(las_fn))
    with open(las_fn, 'wb') as stream:
        write_empty_las(stream)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Write an empty .las file.'
    epilog = 'The point is to pretend that we ran daligner and found no overlaps.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        'las_fn',
        help='Output. A single las file, empty. (12 bytes, actually.)')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
