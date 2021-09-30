module dparasail.log;
import std.stdio;
import std.array : join;

/// Log levels.
enum LogLevel // @suppress(dscanner.style.phobos_naming_convention)
{
    LOG_OFF = 0, ///< All logging disabled.
    LOG_ERROR = 1, ///< Logging of errors only.
    LOG_WARNING = 3, ///< Logging of errors and warnings.
    LOG_INFO = 4, ///< Logging of errors, warnings, and normal but significant events.
    LOG_DEBUG = 5, ///< Logging of all except the most detailed debug events.
    LOG_TRACE = 6 ///< All logging enabled.
}

int verbose;

/// Sets the selected log level.
void setLogLevel(LogLevel level)
{
    verbose = cast(int)level;
}

/// Gets the selected log level.
LogLevel getLogLevel()
{
    return cast(LogLevel) verbose;
}

void log(LogLevel severity, const(char)[] context, string[] contents ...)
{
    final switch(severity)
    {
        case LogLevel.LOG_OFF:
            return;
        case LogLevel.LOG_ERROR:
            stderr.writefln("[E::"~context~"] %s", contents.join);
            return;
        case LogLevel.LOG_WARNING:
            stderr.writefln("[W::"~context~"] %s", contents.join);
            return;
        case LogLevel.LOG_INFO:
            stderr.writefln("[I::"~context~"] %s", contents.join);
            return;
        case LogLevel.LOG_DEBUG:
            stderr.writefln("[D::"~context~"] %s", contents.join);
            return;
        case LogLevel.LOG_TRACE:
            stderr.writefln("[T::"~context~"] %s", contents.join);
            return;
    }
}


/**! Logs an event with severity HTS_LOG_ERROR and default context. Parameters: format, ... */
//#define hts_log_error(...) hts_log(HTS_LOG_ERROR, __func__, __VA_ARGS__)
void logError(const(char)[] ctx, string msg)
{
    string open_error_color = "\x1b[0;31m";
    string close_color      = "\x1b[0m";
    log(LogLevel.LOG_ERROR, ctx, open_error_color, msg, close_color);
}
/**! Logs an event with severity HTS_LOG_WARNING and default context. Parameters: format, ... */
//#define hts_log_warning(...) hts_log(HTS_LOG_WARNING, __func__, __VA_ARGS__)
void logWarning(const(char)[] ctx, string msg)
{
    string open_warning_color = "\x1b[0;33m";
    string close_color        = "\x1b[0m";
    log(LogLevel.LOG_WARNING, ctx, open_warning_color, msg, close_color);
}

/**! Logs an event with severity HTS_LOG_INFO and default context. Parameters: format, ... */
//#define hts_log_info(...) hts_log(HTS_LOG_INFO, __func__, __VA_ARGS__)
void logInfo(const(char)[] ctx, string msg)
{
    string open_info_color = "\x1b[0;32m";
    string close_color     = "\x1b[0m";
    log(LogLevel.LOG_INFO, ctx, open_info_color, msg, close_color);
}

/**! Logs an event with severity HTS_LOG_DEBUG and default context. Parameters: format, ... */
//#define hts_log_debug(...) hts_log(HTS_LOG_DEBUG, __func__, __VA_ARGS__)
void logDebug(const(char)[] ctx, string msg)
{
    string open_debug_color = "\x1b[0;36m";
    string close_color     = "\x1b[0m";
    log(LogLevel.LOG_DEBUG, ctx, open_debug_color, msg, close_color);
}

/**! Logs an event with severity HTS_LOG_TRACE and default context. Parameters: format, ... */
//#define hts_log_trace(...) hts_log(HTS_LOG_TRACE, __func__, __VA_ARGS__)
void logTrace(const(char)[] ctx, string msg)
{
    string open_trace_color = "\x1b[1;36m";
    string close_color     = "\x1b[0m";
    log(LogLevel.LOG_TRACE, ctx, open_trace_color, msg, close_color);
}

///
debug(dhtslib_unittest) unittest
{
    setLogLevel(LogLevel.LOG_TRACE);

    logTrace(__FUNCTION__, "Test: trace");
    logDebug(__FUNCTION__, "Test: debug");
    logInfo(__FUNCTION__,  "Test: info");
    logWarning(__FUNCTION__,"Test: warning");
    logError(__FUNCTION__, "Test: error");
}
