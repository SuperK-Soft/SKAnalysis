# LoadFileList

LoadFileList

## Data

## Configuration

LoadFileList is a general way to specify the input files of ToolChain.
It can accept, in order of precedence:

* `inputFile`, a single filename. This should include the filepath - inputDirectory will NOT be prepended.
* `fileList`, a single filename of a text file containing a list of files to process. Each file should be on a new line. Lines beginning with `#` will be ignored.
* `filePattern`, a pattern by which to match filenames. This may be a glob pattern or a regular expression. If this is to be used, inputDirectory *must* also be specified.
* `inputDirectory`, the path to a directory to search for files matching the filePattern.

Other configuration variables include:
* `FileListName`, this tool will output a vector of strings of filepaths, which will be placed into the CStore. This variable specifies the name with which to retrieve that list. Default is `InputFileList`.
* `useRegex`, when using `filePattern`, whether this represents a regex or a glob pattern.
* `verbosity`, how verbose to be during execution.
