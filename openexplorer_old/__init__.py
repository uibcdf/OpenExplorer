"""
OpenExplorer
Open source Python library to explore the conformation free energy landscape of molecular systems.
"""

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

__documentation_web__ = 'https://www.uibcdf.org/OpenExplorer'
__github_web__ = 'https://github.com/uibcdf/OpenExplorer'
__github_issues_web__ = __github_web__ + '/issues'

from ._pyunitwizard import puw as puw

from .explorer import Explorer
from . import tools
from . import reporters
from .pes import PES
from .ktn import KTN
from . import exploration_campaign
