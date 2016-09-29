import os.path
import re
import shutil
import subprocess
import sys

def get_command_version_string(cmd, regexp, *, prefix="", suffix="", use_stderr=None,
                               encoding=sys.getdefaultencoding(), raise_on_error=False):
    '''Get the version string from a command.

    Arguments:

    cmd: a string or list of strings to run a command that will
    print a version number, typically something like 'mycommand
    --version'.

    regexp: A regular expression including a named capture group with
    a name of 'version' that captures just the version string. For
    example, if the command returns 'mycommand v1.2.3', the regexp
    might be 'v(?P<version>(\\d+\\.)*\\d+)', which would match the
    string '1.2.3' in the 'version' named capture group. If the regexp
    doesn't match the output of the command, or is missing the named
    capture group, an exception is raised.

    Keyword-only arguments:

    prefix, suffix: Prepended/appended to the version string before
    returning.

    use_stderr: If the command is known to print its version to
    standard error instead of standard output, set this to True. If it
    is known to print its version to standard output, set this to
    False. If it is unknown, leave it as None, and the concatenation
    of both (standard output first) will be searched for the version
    string.

    raise_on_error: If False (the default), will return None if an
    error is encountered, including failing to find the command or
    failing to match the regular expression. If True, the exception
    will be raised as normal.

    encoding: Which text encoding to use to read the output of the
    command. Use the system default if not specified.

    '''
    try:
        use_shell = isinstance(cmd, str)
        p = subprocess.Popen(cmd, shell=use_shell, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        if use_stderr is None:
            output = stdout + stderr
        elif use_stderr:
            output = stderr
        else:
            output = stdout
        output = output.decode(encoding)
        m = re.search(regexp, output)
        if m is None:
            raise ValueError("Regular expression did not match command output")
        return prefix + m.group("version") + suffix
    except Exception as ex:
        if raise_on_error:
            raise ex
        else:
            return None

ascp_path = shutil.which("ascp") or os.path.expanduser("~/.aspera/connect/bin/ascp")

# Determine the versions of various programs used
ASCP_VERSION = get_command_version_string([ascp_path, '--version'], 'ascp version\\s+(?P<version>\\S+)', prefix='ascp ')
BOWTIE2_VERSION = get_command_version_string('bowtie2 --version', 'version\\s+(?P<version>\\S+)', prefix='bowtie2 ')
CUFFLINKS_VERSION = get_command_version_string('cufflinks --help', 'cufflinks v(?P<version>\S+)', prefix='cufflinks ')
EPIC_VERSION = get_command_version_string('epic --version', 'epic\\s+(?P<version>\\S+)', prefix='epic ')
FASTQ_TOOLS_VERSION = get_command_version_string('fastq-sort --version', '(?P<version>\\d+(\\.\\d+)*)', prefix='fastq-tools ')
HISAT2_VERSION = get_command_version_string('hisat2 --version', 'version\\s+(?P<version>\\S+)', prefix='hisat2 ')
IDR_VERSION = get_command_version_string('idr --version', '(?P<version>\\d+(\\.\\d+)*)', prefix='IDR ')
KALLISTO_VERSION = get_command_version_string('kallisto', '^kallisto\\s+(?P<version>\\S+)', prefix='kallisto ')
MACS_VERSION = get_command_version_string('macs2 --version', 'macs2\\s+(?P<version>\\S+)', prefix='macs2 ')
SALMON_VERSION = get_command_version_string('salmon --version', 'version\\s+:\\s+(?P<version>\\S+)', prefix='salmon ')
SAMTOOLS_VERSION = get_command_version_string('samtools', 'Version:\\s+(?P<version>\\S+)', prefix='samtools ')
SRATOOLKIT_VERSION = get_command_version_string('fastq-dump --version', ':\\s+(?P<version>\\S+)', prefix='sratoolkit ')
STAR_VERSION = get_command_version_string('STAR --version', 'STAR_(?P<version>\\S+)', prefix='STAR ')
TOPHAT2_VERSION = get_command_version_string('tophat --version', 'TopHat v(?P<version>\\S+)', prefix='tophat ')

# R, BioC, & packages
try:
    from rpy2.robjects import r
    from rpy2.rinterface import RRuntimeError
    R_VERSION = ''.join(r('R.version$version.string'))
except RRuntimeError:
    R_VERSION = None

try:
    from rpy2.robjects import r
    from rpy2.rinterface import RRuntimeError
    BIOC_VERSION = 'Bioconductor ' + "".join(r('as.character(BiocInstaller::biocVersion())'))
except RRuntimeError:
    BIOC_VERSION = None

def R_package_version(pkgname):
    try:
        from rpy2.robjects import r
        from rpy2.rinterface import RRuntimeError
        pkgversion = r('installed.packages()[,"Version"]').rx(pkgname)[0]
        if r['is.na'](pkgversion)[0]:
            raise ValueError("Could not determine package version for {!r}. Maybe the package is not installed?".format(pkgname))
        return ' '.join([pkgname, pkgversion])
    except Exception:
        return None
