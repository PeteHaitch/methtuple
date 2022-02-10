__version_info__ = ('1','7','0')
__version__ = '.'.join(__version_info__)

from . import funcs
from . import mtuple

__all__ = [
    '__version_info__',
    '__version__',
    'funcs',
    'mtuple'
]
