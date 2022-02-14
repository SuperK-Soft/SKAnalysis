// TODO move these to the ToolFrameworkCore Logging class
enum Verbosity
{
    pNONE,    ///< NTag prints no message.
    pERROR,   ///< NTag prints error messages only.
    pWARNING, ///< NTag prints warning and error messages.
    pDEFAULT, ///< NTag prints all messages except debug messages.
    pDEBUG    ///< NTag prints all messages including debug messages.
};
