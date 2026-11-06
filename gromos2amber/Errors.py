
class GromosFormatError(Exception):
    pass

class IllegalArgumentError(ValueError):
    pass

class UnsupportedTopologyError(GromosFormatError):
    """The topology uses a GROMOS feature the converter cannot express.

    A subclass of GromosFormatError so callers that already report format
    errors to the user (the website's system-topology converter) report it.
    """
    pass
