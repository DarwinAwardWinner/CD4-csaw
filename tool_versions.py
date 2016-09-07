import os.path
import re
import shutil
import subprocess
import sys

def get_command_version_string(cmd, regexp, use_stderr=None,
                               encoding=sys.getdefaultencoding()):
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

    use_stderr: If the command is known to print its version to
    standard error instead of standard output, set this to True. If it
    is known to print its version to standard output, set this to
    False. If it is unknown, leave it as None, and the concatenation
    of both (standard output first) will be searched for the version
    string.

    '''
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
    return m.group("version")

# Determine the versions of various programs used
try:
    CUFFLINKS_VERSION = 'cufflinks ' + get_command_version_string('cufflinks --help', 'cufflinks v(?P<version>\S+)')
except Exception:
    CUFFLINKS_VERSION = None

try:
    SAMTOOLS_VERSION = 'samtools ' + get_command_version_string('samtools', 'Version:\\s+(?P<version>\\S+)')
except Exception:
    SAMTOOLS_VERSION = None

try:
    BOWTIE1_VERSION = 'bowtie ' + get_command_version_string('bowtie --version', 'version\\s+(?P<version>\\S+)')
except Exception:
    BOWTIE1_VERSION = None

try:
    BOWTIE2_VERSION = 'bowtie2 ' + get_command_version_string('bowtie2 --version', 'version\\s+(?P<version>\\S+)')
except Exception:
    BOWTIE2_VERSION = None

try:
    TOPHAT2_VERSION = 'tophat ' + get_command_version_string('tophat --version', 'TopHat v(?P<version>\\S+)')
except Exception:
    TOPHAT2_VERSION = None

try:
    HISAT2_VERSION = 'hisat2 ' + get_command_version_string('hisat2 --version', 'version\\s+(?P<version>\\S+)')
except Exception:
    HISAT2_VERSION = None

try:
    BWA_VERSION = 'bwa ' + get_command_version_string('bwa', 'Version:\\s+(?P<version>\\S+)')
except Exception:
    BWA_VERSION = None

try:
    BBMAP_VERSION = 'bbmap ' + get_command_version_string([BBMAP, '--version'], 'BBMap version (?P<version>\\S+)')
except Exception:
    BBMAP_VERSION = None

try:
    STAR_VERSION = 'STAR ' + get_command_version_string('STAR --version', 'STAR_(?P<version>\\S+)')
except Exception:
    STAR_VERSION = None

try:
    SALMON_VERSION = 'salmon ' + get_command_version_string('salmon --version', 'version\\s+:\\s+(?P<version>\\S+)')
except Exception:
    SALMON_VERSION = None

try:
    KALLISTO_VERSION = 'kallisto ' + get_command_version_string('kallisto', '^kallisto\\s+(?P<version>\\S+)')
except Exception:
    KALLISTO_VERSION = None

try:
    ascp_path = shutil.which("ascp") or os.path.expanduser("~/.aspera/connect/bin/ascp")
    ASCP_VERSION = 'ascp ' + get_command_version_string([ascp_path, '--version'], 'ascp version\\s+(?P<version>\\S+)')
except Exception:
    ASCP_VERSION = None

try:
    SRATOOLKIT_VERSION = 'sratoolkit ' + get_command_version_string('fastq-dump --version', ':\\s+(?P<version>\\S+)')
except Exception:
    SRATOOLKIT_VERSION = None

try:
    MACS_VERSION = 'macs2 ' + get_command_version_string('macs2 --version', 'macs2\\s+(?P<version>\\S+)')
except Exception:
    MACS_VERSION = None

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
