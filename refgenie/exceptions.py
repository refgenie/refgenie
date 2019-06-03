__all__ = ["RefgenieError", "MissingGenomeConfigError"]

from refgenconf import CFG_ENV_VARS


class RefgenieError(object):
    pass


class MissingGenomeConfigError(RefgenieError):

    def __init__(self, conf_file=None):
        """
        Create the error message, using optionally an attempt filepath.

        :param str conf_file: path attempted to be used as genome config file
        """
        msg = "Try setting an environment variable to a local config file: " \
              "{}".format(CFG_ENV_VARS)
        if conf_file:
            msg = "Not a file {} -- {}.".format(conf_file, msg)
        super(MissingGenomeConfigError, self).__init__(msg)
