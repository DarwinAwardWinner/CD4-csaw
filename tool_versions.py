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
ASCP_VERSION = get_command_version_string(prefix='ascp ', [ascp_path, '--version'], 'ascp version\\s+(?P<version>\\S+)')
BBMAP_VERSION = get_command_version_string(prefix='bbmap ', [BBMAP, '--version'], 'BBMap version (?P<version>\\S+)')
BOWTIE1_VERSION = get_command_version_string(prefix='bowtie ', 'bowtie --version', 'version\\s+(?P<version>\\S+)')
BOWTIE2_VERSION = get_command_version_string(prefix='bowtie2 ', 'bowtie2 --version', 'version\\s+(?P<version>\\S+)')
BWA_VERSION = get_command_version_string(prefix='bwa ', 'bwa', 'Version:\\s+(?P<version>\\S+)')
CUFFLINKS_VERSION = get_command_version_string(prefix='cufflinks ', 'cufflinks --help', 'cufflinks v(?P<version>\S+)')
EPIC_VERSION = get_command_version_string(prefix='epic ', 'epic --version', 'epic\\s+(?P<version>\\S+)')
FASTQ_TOOLS_VERSION = get_command_version_string(prefix='fastq-tools ', 'fastq-sort --version', '(?P<version>\\d+(\\.\\d+)*)')
HISAT2_VERSION = get_command_version_string(prefix='hisat2 ', 'hisat2 --version', 'version\\s+(?P<version>\\S+)')
IDR_VERSION = get_command_version_string(prefix='IDR ', 'idr --version', '(?P<version>\\d+(\\.\\d+)*)')
KALLISTO_VERSION = get_command_version_string(prefix='kallisto ', 'kallisto', '^kallisto\\s+(?P<version>\\S+)')
MACS_VERSION = get_command_version_string(prefix='macs2 ', 'macs2 --version', 'macs2\\s+(?P<version>\\S+)')
SALMON_VERSION = get_command_version_string(prefix='salmon ', 'salmon --version', 'version\\s+:\\s+(?P<version>\\S+)')
SAMTOOLS_VERSION = get_command_version_string(prefix='samtools ', 'samtools', 'Version:\\s+(?P<version>\\S+)')
SRATOOLKIT_VERSION = get_command_version_string(prefix='sratoolkit ', 'fastq-dump --version', ':\\s+(?P<version>\\S+)')
STAR_VERSION = get_command_version_string(prefix='STAR ', 'STAR --version', 'STAR_(?P<version>\\S+)')
TOPHAT2_VERSION = get_command_version_string(prefix='tophat ', 'tophat --version', 'TopHat v(?P<version>\\S+)')

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
        return ' '.join([pkgname, pkgversion])
    except Exception:
        return None
