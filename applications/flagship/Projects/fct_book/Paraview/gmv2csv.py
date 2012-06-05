#! /usr/bin/env pvpython

"""

Script that converts an ascii or binary GMV file to a CSV file.

Open problem: only convert some arrays, not all.

Usage: gmv2csv.py *.gmv*

Author: Sven Buijssen, sven.buijssen@tu-dortmund.de, 2011

"""

import os;
import sys;
pathname, scriptname = os.path.split(sys.argv[0])
from optparse import OptionParser

try:
  # load python bindings for paraview

  # deprecated
  #from paraview import servermanager;

  from paraview.simple import *;

  # possibly in the future:
  #paraview.simple;

except (ImportError), e:
  # in case the Python bindings are not found, try a few additional paths:
  #path += [ "/sfw/paraview/paraview-3.7.0-Linux-i686/lib/paraview-3.7",
  #          ""
  #        ];
  # or
  #path.append("/sfw/paraview/paraview-3.7.0-Linux-i686/lib/paraview-3.7");
  print
  print scriptname, ": Could not find python bindings for paraview."
  print scriptname, ": Suggested remedy: Add a directory like <..../lib/paraview-3.7>"
  print scriptname, ": to your environment variables PYTHONPATH and LD_LIBRARY_PATH."
  print
  from paraview.simple import *;


def main():
    #
    # Deal with command line options
    #
    usage = "usage: %prog [options] <file list>"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outputdir", dest="outputdirectory",
                      help="store converted files in OUTPUTDIRECTORY,        default: same directory as input files")
    parser.add_option("-q", "--quiet", dest="quietMode", action="store_true", default = 0,
                      help="suppresses screen messages")
    (options, args) = parser.parse_args()
    if len(args) < 1:  # There should only be the input file list left
        parser.error("incorrect number of arguments")
    #output_dir = "";
    output_dir = os.getcwd();
    if options.outputdirectory:
      output_dir = options.outputdirectory;

    #
    # Determine ParaView version for portability issues
    #
    pxm = servermanager.ProxyManager()
    pvVersionString = str(pxm.GetVersionMajor()) + "." + str(pxm.GetVersionMinor()) + "." + str(pxm.GetVersionPatch())
    if not options.quietMode:
      print "Found ParaView " + pvVersionString + "."


    #
    # Loop over all files given on command line
    #
    for infile in args:
      if not os.path.exists(infile):
        raise IOError, '' . join(["\n", scriptname, ": ", "Input file <", infile, "> does not exist.", "\n"]);

      #
      # construct output file name
      #
      if not output_dir == "":
        outfile = os.path.join(output_dir,
                               os.path.basename(infile).replace('.gmv', '.csv'));
      else:
        outfile = infile.replace('.gmv', '.csv');

      #
      # Ensure outfile is not identical to input file so the former does not overwrite the latter
      #
      if outfile == infile:
        outfile += '.csv';
      #print outfile;


      #
      # Connect to ParaView server manager
      #
      if not servermanager.ActiveConnection:
        connection = servermanager.Connect();
        if not connection:
          raise RuntimeError, '' . join(["\n", scriptname, ": ", "Connection could not be established.", "\n"]);


      #
      # Try to establish connection to GMVReader, if necessary try to load it as plugin
      #
      try:
        reader = servermanager.sources.GMVReader(FileNames=[infile]);
      except (AttributeError), e:
        if not options.quietMode:
          print
          print scriptname + ": " + "ParaView has not been built with compile-time GMV Reader support."
          print scriptname + ": " + "Trying to find and load a GMV Reader plugin."

        # Try loading a ParaView Plugin. But to be portable we cannot
        # give absolute paths like:
        # paraview.servermanager.LoadPlugin("/usr/local/paraview/paraview-3.6.1-Linux-i686/bin/plugins/GMVReaderPlugin/libGMVReaderPlugin.so")
        # LoadPlugin() requires an absolute path, though.
        # But a path relative to some path in 'sys.path' would be better and
        # let the script determine the resulting absolute path to the plugin.
        GMVReaderPluginLoaded = 0
        # Prepend a couple of paths to sys.path:
        # The environment variable PV_PLUGIN_PATH and ~/.config/ParaView/ParaView<version>.
        # I'm not sure why that does not happen by default in PV 3.6 and 3.7.
        sys.path.insert(
          0, os.path.join(os.environ['HOME'], ".config", "ParaView",
                          "ParaView" + str(pxm.GetVersionMajor()) + "." + str(pxm.GetVersionMinor()),
                          "Plugins"));
        if os.environ.has_key('PV_PLUGIN_PATH'):
          sep = ';' if sys.platform == 'win32' else ':'
          for path in os.environ['PV_PLUGIN_PATH'].split(sep):
            sys.path.insert(0, path)

        for dir in sys.path:
          # Try <path>/libGMVReaderPlugin.so, <path>/plugins/libGMVReaderPlugin.so.
          for plugindir in [ dir, os.path.join(dir, "plugins") ]:
            pathToGMVReaderPlugin = os.path.join(plugindir, "libGMVReaderPlugin.so")
            print "Looking for GMVReader plugin in " + pathToGMVReaderPlugin + "..."
            if os.path.exists(pathToGMVReaderPlugin):
              paraview.servermanager.LoadPlugin(pathToGMVReaderPlugin, globals())
              GMVReaderPluginLoaded = 1
              if not options.quietMode:
                print scriptname + ": " + "Found GMVReader plugin." + "\n"
              break
          if GMVReaderPluginLoaded:
            break

        if not GMVReaderPluginLoaded:
          raise RuntimeError, '' . join(["\n\n", scriptname, ": ", "No GMVReader Plugin found.", "\n"]);


      #
      # Read file and convert it
      #
      if options.quietMode:
        servermanager.ToggleProgressPrinting()

      reader = servermanager.sources.GMVReader(FileNames=[infile]);
      #   reader.PointArrayStatus = ['velocity', 'u1_x', 'u1_y', 'u2_x', 'u2_y', 'p', 'eps_xx', 'eps_yy', 'eps_xy', 'div_u']
      #   reader.CellArrayStatus = ['material id', 'sigma_xx', 'sigma_yy', 'sigma_xy', 'sigma_zz', 'p_el', 'mises']
      #   reader.PointArrayStatus = ['velocity']
      #   reader.CellArrayStatus = []

      # Does not catch writer errors, sadly.
      #try:
      #  writer.UpdatePipeline();
      #except:
      #  raise RuntimeError, '' . join(["\n", scriptname, ": ", "Conversion of <", infile, "> failed!", "\n"]);

      # Does not catch writer errors either:
      # Hackish solution: If an error occurs, the except branch is *not*
      # executed. So delete output file first (if it exists) and test whether
      # it exists after having called the writer. If it does, print a success
      # message.
      # Does not work as DataSetCSVWriter merges in a number somewhere in the output file.
      # Probably because of the WriteAllTimeSteps option.
      #
      #if os.path.exists(outfile):
      #  os.remove(outfile)
      #
      #writer = servermanager.writers.DataSetCSVWriter(FileName=outfile, Input=reader, WriteAllTimeSteps=0);
      #writer.UpdatePipeline();
      #
      #if os.path.exists(outfile):
      #  if not options.quietMode:
      #    print '' . join([scriptname, ": ", "File <", infile, "> successfully converted to <", outfile, ">."]);

      writer = servermanager.writers.DataSetCSVWriter(FileName=outfile, Input=reader, WriteAllTimeSteps=0);
      writer.UpdatePipeline();


if __name__ == "__main__":
    main()
