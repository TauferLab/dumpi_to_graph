#ifndef D2G_JSONL_FORMAT_H
#define D2G_JSONL_FORMAT_H

// Helpers for writing the JSON Lines event graph partitions (see
// src/EventGraphJsonl.cpp). Kept in a header so test/jsonl_format_test.cpp can
// check them against igraph's GraphML writer.

#include <cmath>
#include <cstdio>
#include <string>

namespace d2g_jsonl {

// JSON string literal. Stops at the first NUL like a C string would, escapes
// what JSON requires, and leaves every other byte as is.
inline void append_json_string( std::string& out, const std::string& value )
{
  out.push_back( '"' );
  for ( char ch : value ) {
    unsigned char c = static_cast<unsigned char>( ch );
    if ( c == 0 ) {
      break;
    }
    switch ( c ) {
      case '"':  out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\b': out += "\\b";  break;
      case '\f': out += "\\f";  break;
      case '\n': out += "\\n";  break;
      case '\r': out += "\\r";  break;
      case '\t': out += "\\t";  break;
      default:
        if ( c < 0x20 ) {
          char buf[8];
          std::snprintf( buf, sizeof(buf), "\\u%04x", c );
          out += buf;
        } else {
          out.push_back( ch );
        }
    }
  }
  out.push_back( '"' );
}

// Shortest-safe round-trip representation of a double
inline void append_json_double( std::string& out, double value )
{
  if ( std::isnan( value ) ) {
    out += "\"nan\"";
  } else if ( std::isinf( value ) ) {
    out += value > 0 ? "\"inf\"" : "\"-inf\"";
  } else {
    char buf[40];
    std::snprintf( buf, sizeof(buf), "%.17g", value );
    out += buf;
  }
}

inline void append_uint( std::string& out, unsigned long long value )
{
  char buf[32];
  std::snprintf( buf, sizeof(buf), "%llu", value );
  out += buf;
}

inline void append_int( std::string& out, long long value )
{
  char buf[32];
  std::snprintf( buf, sizeof(buf), "%lld", value );
  out += buf;
}

} // namespace d2g_jsonl

#endif // D2G_JSONL_FORMAT_H
