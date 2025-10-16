#include "gdal.h"
#include "gdal_priv.h"
#include "nlohmann/json.hpp"

// C API
#ifdef HAS_GEOS
#include "geos_c.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <algorithm>

using json = nlohmann::json;
using namespace std;

class RasterProcessor {
private:
    GDALDataset* dataset;
    int width, height;
    vector<uint8_t> maskData;
    double geoTransform[6];
    bool hasGeoTransform;
    
public:
    RasterProcessor() : dataset(nullptr), hasGeoTransform(false) {
        GDALAllRegister();
    }
    
    ~RasterProcessor() {
        if (dataset) {
            GDALClose(dataset);
        }
    }
    
    bool loadRaster(const string& filename) {
        dataset = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
        if (!dataset) {
            cerr << "Не удалось открыть файл: " << filename << endl;
            return false;
        }
        
        width = dataset->GetRasterXSize();
        height = dataset->GetRasterYSize();
        hasGeoTransform = (dataset->GetGeoTransform(geoTransform) == CE_None);
        
        cout << "Загружен: " << filename << " (" << width << "x" << height << ")" << endl;
        
        return loadMaskData();
    }
    
    #ifdef HAS_GEOS
    GEOSGeometry* getValidGeometry() {
        if (maskData.empty()) {
            return nullptr;
        }
        
        cout << "Создание геометрии из маски..." << endl;
        
        // Границы непрозрачных данных
        int dataMinX = width, dataMinY = height, dataMaxX = 0, dataMaxY = 0;
        bool foundData = false;
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                if (isOpaque(x, y)) {
                    foundData = true;
                    dataMinX = min(dataMinX, x);
                    dataMinY = min(dataMinY, y);
                    dataMaxX = max(dataMaxX, x);
                    dataMaxY = max(dataMaxY, y);
                }
            }
        }
        
        if (!foundData) {
            cout << "Непрозрачные данные не найдены" << endl;
            return nullptr;
        }
        
        cout << "Границы данных: [" << dataMinX << "," << dataMinY << "] - [" 
             << dataMaxX << "," << dataMaxY << "]" << endl;
        
        // Преобразуем в географические координаты
        double ulx, uly, lrx, lry;
        pixelToGeo(dataMinX, dataMinY, ulx, uly);
        pixelToGeo(dataMaxX + 1, dataMaxY + 1, lrx, lry);
        
        cout << "Географические границы: [" << ulx << "," << uly << "] - [" 
             << lrx << "," << lry << "]" << endl;
        
        // Создаем полигон через GEOS C API
        GEOSCoordSequence* coordSeq = GEOSCoordSeq_create(5, 2);
        
        GEOSCoordSeq_setX(coordSeq, 0, ulx);
        GEOSCoordSeq_setY(coordSeq, 0, uly);
        
        GEOSCoordSeq_setX(coordSeq, 1, lrx);
        GEOSCoordSeq_setY(coordSeq, 1, uly);
        
        GEOSCoordSeq_setX(coordSeq, 2, lrx);
        GEOSCoordSeq_setY(coordSeq, 2, lry);
        
        GEOSCoordSeq_setX(coordSeq, 3, ulx);
        GEOSCoordSeq_setY(coordSeq, 3, lry);
        
        GEOSCoordSeq_setX(coordSeq, 4, ulx);
        GEOSCoordSeq_setY(coordSeq, 4, uly);
        
        GEOSGeometry* ring = GEOSGeom_createLinearRing(coordSeq);
        GEOSGeometry* polygon = GEOSGeom_createPolygon(ring, nullptr, 0);
        
        return polygon;
    }
    #endif
    
    void printDetailedInfo() {
        if (!dataset) return;
        
        cout << "\nДетальная информация о растре:" << endl;
        cout << "Размер: " << width << "x" << height << endl;
        cout << "Каналы: " << dataset->GetRasterCount() << endl;
        
        for (int i = 1; i <= dataset->GetRasterCount(); ++i) {
            GDALRasterBand* band = dataset->GetRasterBand(i);
            GDALColorInterp colorType = band->GetColorInterpretation();
            const char* colorName = GDALGetColorInterpretationName(colorType);
            cout << "  Канал " << i << ": " << colorName << endl;
        }
        
        // Статистика по маске
        int opaqueCount = 0;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                if (isOpaque(x, y)) {
                    opaqueCount++;
                }
            }
        }
        
        cout << "Непрозрачных пикселей: " << opaqueCount 
             << " (" << (opaqueCount * 100.0 / (width * height)) << "%)" << endl;
        
        // Информация о геотрансформации
        if (hasGeoTransform) {
            cout << "Геотрансформация: [" << geoTransform[0] << ", " << geoTransform[1] 
                 << ", " << geoTransform[2] << ", " << geoTransform[3] 
                 << ", " << geoTransform[4] << ", " << geoTransform[5] << "]" << endl;
        } else {
            cout << "Геотрансформация не найдена" << endl;
        }
    }
    
private:
    bool loadMaskData() {
        int alphaBand = findAlphaBand();
        if (alphaBand == -1) {
            cerr << "Альфа-канал не найдена" << endl;
            return false;
        }
        
        cout << "Загрузка альфа-канала (канал " << alphaBand << ")..." << endl;
        
        GDALRasterBand* band = dataset->GetRasterBand(alphaBand);
        maskData.resize(width * height);
        
        CPLErr err = band->RasterIO(GF_Read, 0, 0, width, height,
                                  maskData.data(), width, height,
                                  GDT_Byte, 0, 0);
        
        return err == CE_None;
    }
    
    int findAlphaBand() {
        int bandCount = dataset->GetRasterCount();
        
        for (int i = 1; i <= bandCount; ++i) {
            GDALRasterBand* band = dataset->GetRasterBand(i);
            GDALColorInterp colorInterp = band->GetColorInterpretation();
            
            if (colorInterp == GCI_AlphaBand) {
                return i;
            }
        }
        
        if (bandCount >= 4) {
            return bandCount;
        }
        
        return -1;
    }
    
    bool isOpaque(int x, int y) {
        return maskData[y * width + x] == 255;
    }
    
    void pixelToGeo(int x, int y, double& geoX, double& geoY) {
        if (hasGeoTransform) {
            geoX = geoTransform[0] + x * geoTransform[1] + y * geoTransform[2];
            geoY = geoTransform[3] + x * geoTransform[4] + y * geoTransform[5];
        } else {
            geoX = x;
            geoY = height - y;
        }
    }
};

#ifdef HAS_GEOS
string geometryToGeoJSON(GEOSGeometry* geometry) {
    if (!geometry) {
        return R"({"type": "FeatureCollection", "features": []})";
    }
    
    // Получаем bounding box геометрии
    double minX = 0.0, minY = 0.0, maxX = 100.0, maxY = 100.0;
    
    // Простой способ - используем envelope для получения bounding box
    GEOSGeometry* envelope = GEOSEnvelope(geometry);
    if (envelope) {
        // Получаем внешнее кольцо envelope
        const GEOSGeometry* ring = GEOSGetExteriorRing(envelope);
        if (ring) {
            const GEOSCoordSequence* coordSeq = GEOSGeom_getCoordSeq(ring);
            if (coordSeq) {
                unsigned int size;
                if (GEOSCoordSeq_getSize(coordSeq, &size)) {
                    // Берем первую точку (min) и третью точку (max) для bounding box
                    GEOSCoordSeq_getX(coordSeq, 0, &minX);
                    GEOSCoordSeq_getY(coordSeq, 0, &minY);
                    GEOSCoordSeq_getX(coordSeq, 2, &maxX);
                    GEOSCoordSeq_getY(coordSeq, 2, &maxY);
                }
            }
        }
        GEOSGeom_destroy(envelope);
    }
    
    json geojson = {
        {"type", "FeatureCollection"},
        {"features", {
            {
                {"type", "Feature"},
                {"properties", {
                    {"name", "Intersection Area"}
                }},
                {"geometry", {
                    {"type", "Polygon"},
                    {"coordinates", {{
                        {minX, minY},
                        {maxX, minY},
                        {maxX, maxY},
                        {minX, maxY},
                        {minX, minY}
                    }}}
                }}
            }
        }}
    };
    
    return geojson.dump(4);
}
#endif

int main() {
    #ifdef HAS_GEOS
    // Инициализируем GEOS
    initGEOS(nullptr, nullptr);
    #endif
    
    try {
        cout << "=== Анализатор пересечения растров (GEOS C API) ===" << endl;
        
        // Обрабатываем первое изображение
        RasterProcessor processor1;
        if (!processor1.loadRaster("orto1.tif")) {
            cerr << "Ошибка загрузки orto1.tif" << endl;
            return 1;
        }
        processor1.printDetailedInfo();
        
        // Обрабатываем второе изображение  
        RasterProcessor processor2;
        if (!processor2.loadRaster("orto2.tif")) {
            cerr << "Ошибка загрузки orto2.tif" << endl;
            return 1;
        }
        processor2.printDetailedInfo();
        
        #ifdef HAS_GEOS
        cout << "\nВычисление пересечения..." << endl;
        
        GEOSGeometry* geometry1 = processor1.getValidGeometry();
        GEOSGeometry* geometry2 = processor2.getValidGeometry();
        
        if (geometry1 && geometry2) {
            cout << "Геометрии созданы успешно" << endl;
            
            // Вычисляем пересечение
            GEOSGeometry* intersection = GEOSIntersection(geometry1, geometry2);
            
            if (intersection && !GEOSisEmpty(intersection)) {
                cout << "Пересечение найдено!" << endl;
                
                string geojson = geometryToGeoJSON(intersection);
                ofstream file("intersection_obchaja_2.geojson");
                file << geojson;
                file.close();
                
                cout << "Файл intersection_obchaja_2.geojson создан!" << endl;
                
                GEOSGeom_destroy(intersection);
            } else {
                cout << "Пересечение не найдено или пустое" << endl;
                
                json emptyGeojson = {
                    {"type", "FeatureCollection"},
                    {"features", json::array()}
                };
                
                ofstream file("intersection_obchaja_2.geojson");
                file << emptyGeojson.dump(4);
                file.close();
            }
            
            GEOSGeom_destroy(geometry1);
            GEOSGeom_destroy(geometry2);
        } else {
            cerr << "Не удалось создать геометрии" << endl;
            if (!geometry1) cerr << "  - Геометрия 1 не создана" << endl;
            if (!geometry2) cerr << "  - Геометрия 2 не создана" << endl;
        }
        #else
        cout << "GEOS не доступен" << endl;
        #endif
        
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    
    #ifdef HAS_GEOS
    // Завершаем GEOS
    finishGEOS();
    #endif
    
    cout << "\nПрограмма завершена" << endl;
    return 0;

}
