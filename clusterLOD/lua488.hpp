// Termm-Fall 2024

#ifndef LUA_HPP
#define LUA_HPP

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32  // Windows 平台
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#else          // 非 Windows 平台
#include <lua-5.4.6/src/lua.h>
#include <lua-5.4.6/src/lualib.h>
#include <lua-5.4.6/src/lauxlib.h>
#endif

#ifdef __cplusplus
}
#endif

#endif
