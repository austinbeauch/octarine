"""Validation application. Runs locally at http://localhost:9909/app"""
import sys
import daomop.viewer


if __name__ == "__main__":

    from optparse import OptionParser

    usage = "usage: %prog [options] cmd [args]"
    optprs = OptionParser(usage=usage, version='%%prog')
    optprs.add_option("--debug", dest="debug", default=False, action="store_true",
                      help="Enter the pdb debugger on main()")
    optprs.add_option("--host", dest="host", metavar="HOST",
                      default='localhost',
                      help="Listen on HOST for connections")
    optprs.add_option("--log", dest="logfile", metavar="FILE",
                      help="Write logging output to FILE")
    optprs.add_option("--loglevel", dest="loglevel", metavar="LEVEL",
                      type='int', default=daomop.viewer.logging.INFO,
                      help="Set logging level to LEVEL")
    optprs.add_option("--port", dest="port", metavar="PORT",
                      type=int, default=9909,
                      help="Listen on PORT for connections")
    optprs.add_option("--stderr", dest="logstderr", default=False,
                      action="store_true",
                      help="Copy logging also to stderr")

    (options, args) = optprs.parse_args(sys.argv[1:])

    if options.debug:
        import pdb
        pdb.run('daomop.viewer.main(options)')

    daomop.viewer.main(options)
