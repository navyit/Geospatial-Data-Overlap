cmake_minimum_required(VERSION 3.12)
project(RasterIntersection)

set(CMAKE_CXX_STANDARD 17)

# Явно указываем пути к vcpkg
set(VCPKG_ROOT "C:/vcpkg/installed/x64-windows")
include_directories("${VCPKG_ROOT}/include")

find_package(GDAL REQUIRED)
find_package(nlohmann_json REQUIRED)

# Ищем GEOS C library
find_library(GEOS_C_LIBRARY NAMES geos_c
    PATHS "${VCPKG_ROOT}/lib")

add_executable(raster_intersection main2.cpp)

# Подключаем библиотеки
target_link_libraries(raster_intersection 
    GDAL::GDAL 
    nlohmann_json::nlohmann_json
)

# Добавляем include директории
target_include_directories(raster_intersection PRIVATE 
    "${VCPKG_ROOT}/include"
)

# Добавляем GEOS C если найден
if(GEOS_C_LIBRARY)
    target_link_libraries(raster_intersection ${GEOS_C_LIBRARY})
    target_compile_definitions(raster_intersection PRIVATE HAS_GEOS=1)
    message(STATUS "✅ GEOS C подключен: ${GEOS_C_LIBRARY}")
else()
    message(WARNING "❌ GEOS C не найден")
endif()

message(STATUS "✅ Конфигурация завершена!")
