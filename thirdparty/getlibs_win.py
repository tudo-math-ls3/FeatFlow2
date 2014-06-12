import os
import sys
import subprocess as sp

# third-party package dictionary
# The first part of the value is the download url.
# The second part is the name of a file that is checked to determine whether a
# particular package has already been unpacked.
# Note: The keys of the dictionary coincide with the unpack directories.
pack = {} # [url, testfile]
pack["BLAS"] = [
  "http://www.netlib.org/blas/blas.tgz",
  "xerbla.f"
]
pack["lapack-3.5.0"] = [
  "http://www.netlib.org/lapack/lapack-3.5.0.tgz",
  "LICENSE"
]
pack["AMD"] = [
  "http://www.cise.ufl.edu/research/sparse/amd/AMD-2.3.1.tar.gz",
  "README.txt"
]
pack["SuiteSparse_config"] = [
  "http://www.cise.ufl.edu/research/sparse/SuiteSparse_config/SuiteSparse_config-4.2.1.tar.gz",
  "README.txt"
]
pack["UMFPACK"] = [
  "http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK-5.6.2.tar.gz",
  "README.txt"
]

# download a file from the internet
def download(url,filename):
  if sys.version_info[0] < 3:
    # This code works only in Python 2.x; it does not work for Pyton 3.x
    import urllib
    urllib.urlretrieve(url, os.path.join(".",filename))
    return True
  else:
    import urllib.request
    import shutil
    try:
      with urllib.request.urlopen(url) as response:
        with open(os.path.join(".",filename), 'wb') as out_file:
          shutil.copyfileobj(response, out_file)
    except urllib.error.HTTPError as e:
      print ()
      print("ERROR: Failed to download '" + filename + "'")
      print("REASON: " + str(e.code) + " " + e.reason)
      return False
    return True

# unpack a (gzipped) tarball
def unpack(filename, outdir):
  import tarfile
  try:
    archive = tarfile.open(filename, "r")
    archive.extractall(outdir)
  except FileNotFoundError:
    print ()
    print("ERROR: Archive not found: '" + filename + "'")
    return False
  return True

##################################### MAIN #####################################

print("")
print("Welcome to the FeatFlow2 third-party download script")

# check args
purge = ("--purge" in sys.argv[1:])
force_download = ("--force-download" in sys.argv[1:])
force_unpack   = ("--force-unpack" in sys.argv[1:])

# done, fail, skip counters
done_dl = 0
fail_dl = 0
skip_dl = 0
done_up = 0
fail_up = 0
skip_up = 0

# purge first?
if purge:
  sp.call("purgefiles_win.cmd")

# loop over all URLs
for name in pack:
  url = pack[name][0]
  filename = url.split("/")[-1]
  dirpath = os.path.abspath(os.path.join(".", name))
  filepath = os.path.abspath(os.path.join(".",filename))
  testfile = os.path.abspath(os.path.join(".",name,pack[name][1]))

  # download file
  if (force_download) or (not os.path.isfile(filepath)):
    print ("Downloading '" + filename + "'...")
    if not download(url, filepath):
      fail_dl += 1
      continue
    done_dl += 1
  else:
    print ("INFO: Archive '" + filepath + "' already exists; skipping download")
    skip_dl += 1

  # unpack file
  if (force_unpack) or (not os.path.isfile(testfile)):
    print("Unpacking '" + filename + "'...")
    if not unpack(filepath, "."):
      fail_up += 1
      continue
    done_up += 1
  else:
    print ("INFO: Archive '" + filepath + "' already unpacked; skipping unpack")
    skip_up += 1

# print summary
print("")
print("Downloads: " + str(done_dl) + " succeeded, " + str(fail_dl) + " failed, " + str(skip_dl) + " skipped")
print("Unpacks..: " + str(done_up) + " succeeded, " + str(fail_up) + " failed, " + str(skip_up) + " skipped")
