#include <QCoreApplication>
#include <QDebug>
#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "ogr_core.h"
#include "ogr_api.h"

const double PI = 3.14159265;
const double ANGLE_OFFSET[5] = {-90.0, -45.0, 0.0, 45.0, 90.0};

int run(const char *inBrat, const char *inDem, const char *outPoints, const char *outPonds);
double addDegrees(double base, double addValue);
double calcAzimuth(double startX, double startY, double endX, double endY);
int calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY);
double getRasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord);
int loadRasterValuesFromPoint(const char *rasterPath, const char *pointPath, int rows, int cols, double transform[]);

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    const char *bratIn = "E:/etal/Projects/NonLoc/Beaver_Modeling/02_Data/BRAT/TempleFk";
    const char *demIn = "E:/etal/Projects/NonLoc/Beaver_Modeling/02_Data/DEM/NED_10m/TempleFk_WS/templefk_10m_ws.tif";
    const char *pointsOut = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/pondpoints.shp";
    const char *pondsOut = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/pondstorage.tif";

    qDebug()<<"variables declared starting run";

    run(bratIn, demIn, pointsOut, pondsOut);

    return a.exec();
    //return 0;
}

int run(const char *inBrat, const char *inDem, const char *outPoints, const char *outPonds)
{
    qDebug()<<"in run";
    GDALAllRegister();
    OGRRegisterAll();
    qDebug()<<"gdal registered";

    const char *startPointsOut = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/startpoints";
    const char *endPointsOut = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/endpoints";
    const char *polygonsOut = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/polygons";

    GDALDataset *pDem, *pPondDepth;
    OGRDataSource *pBrat, *pPondPoints, *pStartPoints, *pEndPoints, *pPolygons;

    GDALDriver *pDriverTiff;
    pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");
    qDebug()<<"gdal driver created";

    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");
    pStartPoints = pDriverShp->CreateDataSource(startPointsOut, NULL);
    pEndPoints = pDriverShp->CreateDataSource(endPointsOut, NULL);
    pPolygons = pDriverShp->CreateDataSource(polygonsOut, NULL);
    qDebug()<<"ogr driver created";

    qDebug()<<"gdal dataset opened";
    pBrat = pDriverShp->CreateDataSource(inBrat, NULL);
    qDebug()<<"ogr dataset opened";

    qDebug()<<pBrat->GetName()<<pBrat->GetLayerCount()<<pBrat->GetLayer(0)->GetName();

    OGRLayer *pBratLayer = pBrat->GetLayer(0);
    int nFeatures = pBratLayer->GetFeatureCount();
    OGRFeature *pBratFeature, *pStartFeature, *pEndFeature, *pPolygonFeature;
    OGRGeometry *pGeom;

    pStartPoints->CreateLayer("BratStartPoints", pBratLayer->GetSpatialRef(), wkbPoint, NULL);
    pEndPoints->CreateLayer("BratEndPoints", pBratLayer->GetSpatialRef(), wkbPoint, NULL);
    pPolygons->CreateLayer("PondSearchPolygons", pBratLayer->GetSpatialRef(), wkbPolygon, NULL);

    OGRLayer *pStartLayer, *pEndLayer, *pPolyLayer;
    pStartLayer = pStartPoints->GetLayer(0);
    pEndLayer = pEndPoints->GetLayer(0);
    pPolyLayer = pPolygons->GetLayer(0);
    OGRFieldDefn field("dam_elev", OFTReal);
    pPolyLayer->CreateField(&field);
    pStartFeature = OGRFeature::CreateFeature(pStartLayer->GetLayerDefn());
    pEndFeature = OGRFeature::CreateFeature(pEndLayer->GetLayerDefn());
    pPolygonFeature = OGRFeature::CreateFeature(pPolyLayer->GetLayerDefn());
    qDebug()<<"points layer(s) added";

    OGRPoint point1, point2;
    int nPoints;
    double azimuthStart, azimuthCurrent, distance, pointX, pointY;
    QVector<double> xCoords(5), yCoords(5);

    for (int i=0; i<nFeatures; i++)
    {
        qDebug()<<i;
        pBratFeature = pBratLayer->GetFeature(i);
        qDebug()<<"brat layer retrieved";
        pGeom = pBratFeature->GetGeometryRef();
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        nPoints = pBratLine->getNumPoints();
        qDebug()<<"nPoints "<<nPoints;
        qDebug()<<pBratLine->getX(0)<<pBratLine->getY(0);
        point1.setX(pBratLine->getX(0));
        point1.setY(pBratLine->getY(0));
        qDebug()<<point1.getX()<<point1.getY();
        pStartFeature->SetGeometry(&point1);
        qDebug()<<"geom set";
        pStartLayer->CreateFeature(pStartFeature);
        point2.setX(pBratLine->getX(nPoints-1));
        point2.setY(pBratLine->getY(nPoints-1));
        pEndFeature->SetGeometry(&point2);
        pEndLayer->CreateFeature(pEndFeature);
        double damElevation, damHeight;
        damElevation = getRasterValueAtPoint(inDem, point2.getX(), point2.getY());
        damHeight = 1.0;
        qDebug()<<damElevation<<damHeight;
        damElevation += damHeight;
        qDebug()<<"dam elevation = "<<damElevation;
        pPolygonFeature->SetField("dam_elev", damElevation);
        qDebug()<<"field set";

        azimuthStart = calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        distance = 25.0;
        pointX=0, pointY=0;
        qDebug()<<"setting poly and linear ring";
        OGRPolygon pPondPoly;
        OGRLinearRing pRing;
        OGRPoint pPolyPoint;

        qDebug()<<"starting poly loop";
        for (int i=0; i<5; i++)
        {
            azimuthCurrent = addDegrees(azimuthStart, ANGLE_OFFSET[i]);
            calcCoords(point2.getX(), point2.getY(), azimuthCurrent, distance, xCoords[i], yCoords[i]);
            qDebug()<<"coords calced"<<xCoords[i]<<yCoords[i];
            pPolyPoint.setX(xCoords[i]);
            pPolyPoint.setY(yCoords[i]);
            qDebug()<<"point created";
            pRing.addPoint(&pPolyPoint);
            qDebug()<<"point added";
        }
        qDebug()<<"finishing loop";
        pRing.addPoint(xCoords[0], yCoords[0]);
        qDebug()<<"final point added";
        pPondPoly.addRing(&pRing);
        qDebug()<<"ring added";
        pPolygonFeature->SetGeometry(&pPondPoly);
        qDebug()<<"poly geom set";
        pPolyLayer->CreateFeature(pPolygonFeature);
        qDebug()<<"feature added";
    }

    pDem = (GDALDataset*) GDALOpen(inDem, GA_ReadOnly);
    int nRows, nCols;
    nCols = pDem->GetRasterXSize();
    nRows = pDem->GetRasterYSize();
    double transform[6];
    pDem->GetGeoTransform(transform);

    const char *waterRas = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/waterras.tif";
    const char *waterPts = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/pond_pts.shp";

    loadRasterValuesFromPoint(waterRas, waterPts, nRows, nCols, transform);

    OGRFeature::DestroyFeature(pBratFeature);
    OGRFeature::DestroyFeature(pStartFeature);
    OGRFeature::DestroyFeature(pEndFeature);
    OGRFeature::DestroyFeature(pPolygonFeature);

    OGRDataSource::DestroyDataSource(pBrat);
    OGRDataSource::DestroyDataSource(pStartPoints);
    OGRDataSource::DestroyDataSource(pEndPoints);
    OGRDataSource::DestroyDataSource(pPolygons);

    GDALClose(pDem);
    qDebug()<<"done";

    return 0;
}

double addDegrees(double base, double addValue)
{
    double value, remainder;

    value = base + addValue;

    if (value > 360.0)
    {
        remainder = value - 360.0;
        value = remainder;
    }
    else if (value < 0.0)
    {
        remainder = fabs(value);
        value = 360.0 - remainder;
    }

    return value;
}

double calcAzimuth(double startX, double startY, double endX, double endY)
{
    double theta, azimuth;

    if (startX > endX && startY > endY)
    {
        theta = atan2((startY-endY),(startX-endX)) * 180.0 / PI;
        azimuth = 180.0 + theta;
    }
    else if (startX < endX && startY > endY)
    {
        theta = atan2((startY-endY),(endX-startX)) * 180.0 / PI;
        azimuth = 360.0 - theta;
    }
    else if (startX < endX && startY < endY)
    {
        theta = atan2((endY-startY),(endX-startX)) * 180.0 / PI;
        azimuth = theta;
    }
    else if (startX > endX && startY < endY)
    {
        theta = atan2((startX-endX),(endY-startY)) * 180.0 / PI;
        azimuth = 90.0 + theta;
    }
    else
    {
        qDebug()<<"AZIMUTH ERROR: may be a straight line";
    }

    return azimuth;
}

int calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY)
{
    double deltaX, deltaY, theta;

    if (azimuth > 0.0 && azimuth < 90.0)
    {
        theta = azimuth;
        deltaY = sin(theta*(PI/180.0))*distance;
        deltaX = cos(theta*(PI/180.0))*distance;
    }
    else if (azimuth > 90.0 && azimuth < 180.0)
    {
        theta = azimuth - 90.0;
        deltaX = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaY = cos(theta*(PI/180.0))*distance;
    }
    else if (azimuth > 180.0 && azimuth < 270.0)
    {
        theta = azimuth - 180.0;
        deltaY = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaX = cos(theta*(PI/180.0))*distance*(-1.0);
    }
    else if (azimuth > 270.0 && azimuth < 360.0)
    {
        theta = 360.0 - azimuth;
        deltaY = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaX = cos(theta*(PI/180.0))*distance;
    }
    else
    {

    }

    newX = startX + deltaX;
    newY = startY + deltaY;

    return 0;
}

double getRasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord)
{
    GDALDataset *pRaster;
    double transform[6];
    double value, xOffset, yOffset, xDiv, yDiv;
    int row, col;

    pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);

    pRaster->GetGeoTransform(transform);

    xOffset = xCoord - transform[0];
    yOffset = transform[3] - yCoord;

    xDiv = xOffset/transform[1];
    yDiv = yOffset/transform[1];

    row = floor(yDiv);
    col = floor(xDiv);

    float *rasVal = (float*) CPLMalloc(sizeof(float)*1);

    pRaster->GetRasterBand(1)->RasterIO(GF_Read, col, row, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);

    value = *rasVal;

    GDALClose(pRaster);
    CPLFree(rasVal);

    return value;
}

int loadRasterValuesFromPoint(const char *rasterPath, const char *pointPath, int rows, int cols, double transform[6])
{
    GDALDataset *pRaster;
    GDALDriver *pDriver;
    OGRDataSource *pPoints;

    double noData = -9999;

    pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    pRaster = pDriver->Create(rasterPath,cols, rows, 1, GDT_Float32, NULL);
    pRaster->SetGeoTransform(transform);
    pRaster->GetRasterBand(1)->Fill(noData);
    pRaster->GetRasterBand(1)->SetNoDataValue(noData);

    pPoints = OGRSFDriverRegistrar::Open(pointPath);
    OGRLayer *layer = pPoints->GetLayer(0);
    int nFeatures = layer->GetFeatureCount();
    qDebug()<<"pond points"<<nFeatures;

}
