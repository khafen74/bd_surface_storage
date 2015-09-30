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
int createDamPoints(const char *bratPath, const char *pointPath);
int getRasterCol(double transform[6], double xCoord);
int getRasterRow(double transform[6], double yCoord);
double getRasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord);
int loadRasterValuesFromPoint(const char *rasterPath, const char *pointPath, int rows, int cols, double transform[]);

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    const char *bratIn = "E:/etal/Projects/NonLoc/Beaver_Modeling/02_Data/BRAT/TempleFk";
    const char *demIn = "E:/etal/Projects/NonLoc/Beaver_Modeling/02_Data/DEM/LiDaR_1m/fme450000.tif";
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

    const char *pondPointsOut = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/pondpoints";
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
    pPondPoints = pDriverShp->CreateDataSource(pondPointsOut, NULL);
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
    OGRGeometry *pGeom;

    pPondPoints->CreateLayer("PondPoints", pBratLayer->GetSpatialRef(), wkbPoint, NULL);
    pStartPoints->CreateLayer("BratStartPoints", pBratLayer->GetSpatialRef(), wkbPoint, NULL);
    pEndPoints->CreateLayer("BratEndPoints", pBratLayer->GetSpatialRef(), wkbPoint, NULL);
    pPolygons->CreateLayer("PondSearchPolygons", pBratLayer->GetSpatialRef(), wkbPolygon, NULL);

    OGRLayer *pStartLayer, *pEndLayer, *pPolyLayer, *pPondLayer;
    pPondLayer = pPondPoints->GetLayer(0);
    pStartLayer = pStartPoints->GetLayer(0);
    pEndLayer = pEndPoints->GetLayer(0);
    pPolyLayer = pPolygons->GetLayer(0);
    OGRFieldDefn field("dam_elev", OFTReal);
    pPolyLayer->CreateField(&field);
    qDebug()<<"points layer(s) added";

    OGRPoint point1, point2;
    int nPoints, nDamCount;
    double azimuthStart, azimuthCurrent, distance, pointX, pointY, length, damDens, spacing, slope;
    QVector<double> xCoords(5), yCoords(5);

    for (int i=0; i<nFeatures; i++)
    {
        OGRFeature *pBratFeature, *pStartFeature, *pEndFeature, *pPondFeature;
        pPondFeature = OGRFeature::CreateFeature(pPondLayer->GetLayerDefn());
        pStartFeature = OGRFeature::CreateFeature(pStartLayer->GetLayerDefn());
        pEndFeature = OGRFeature::CreateFeature(pEndLayer->GetLayerDefn());

        pBratFeature = pBratLayer->GetFeature(i);

        pGeom = pBratFeature->GetGeometryRef();
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        nPoints = pBratLine->getNumPoints();

        point1.setX(pBratLine->getX(0));
        point1.setY(pBratLine->getY(0));

        pStartFeature->SetGeometry(&point1);

        pStartLayer->CreateFeature(pStartFeature);
        point2.setX(pBratLine->getX(nPoints-1));
        point2.setY(pBratLine->getY(nPoints-1));
        pEndFeature->SetGeometry(&point2);
        pEndLayer->CreateFeature(pEndFeature);


        const char *densName = "e_DamDens";
        const char *slopeName = "iGeo_Slope";
        length = pBratLine->get_Length();
        damDens = pBratFeature->GetFieldAsDouble(densName);
        slope = pBratFeature->GetFieldAsDouble(slopeName);
        nDamCount = 0;
        nDamCount = round(length * (damDens/1000.0));
        if (nDamCount > 0)
        {
            spacing = length/(nDamCount*1.0);
        }
        else
        {
            spacing = 0.0;
        }

        if (damDens <= 0.0 && nDamCount > 0)
        {
            qDebug()<<"error: dam dens "<<damDens<<nDamCount;
        }
        //qDebug()<<"space, length, dens, count"<<spacing<<length<<damDens<<nDamCount;

        double damElevation, damHeight, pointDist;
        damHeight = 1.0;

        azimuthStart = calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        distance = (damHeight/slope)*2.0;
        if (distance > spacing)
        {
            distance = spacing;
        }
        pointX=0, pointY=0;

        for (int j=0; j<nDamCount; j++)
        {
            OGRPoint damPoint;
            OGRFeature *pPolygonFeature;
            pPolygonFeature = OGRFeature::CreateFeature(pPolyLayer->GetLayerDefn());
            pointDist = length - (spacing * (j * 1.0));
            //qDebug()<<j+1<<" of "<<nDamCount<<" point dist "<<pointDist;
            pBratLine->Value(pointDist, &damPoint);
            //qDebug()<<"dam point retrieved";
            damElevation = getRasterValueAtPoint(inDem, damPoint.getX(), damPoint.getY());
            //qDebug()<<"dam elevation retrieved";
            damElevation += damHeight;

            pPolygonFeature->SetField("dam_elev", damElevation);
            //qDebug()<<"field set"<<damElevation;
            OGRPolygon pPondPoly;
            OGRLinearRing pRing;
            OGRPoint pPolyPoint;
            //qDebug()<<"starting poly loop";

            for (int k=0; k<5; k++)
            {
                azimuthCurrent = addDegrees(azimuthStart, ANGLE_OFFSET[k]);
                calcCoords(damPoint.getX(), damPoint.getY(), azimuthCurrent, distance, xCoords[k], yCoords[k]);

                pPolyPoint.setX(xCoords[k]);
                pPolyPoint.setY(yCoords[k]);

                pRing.addPoint(&pPolyPoint);
            }

            pRing.addPoint(xCoords[0], yCoords[0]);
            pPondPoly.addRing(&pRing);
            pPolygonFeature->SetGeometry(&pPondPoly);
            pPolyLayer->CreateFeature(pPolygonFeature);
            OGRFeature::DestroyFeature(pPolygonFeature);
        }

        OGRFeature::DestroyFeature(pStartFeature);
        OGRFeature::DestroyFeature(pEndFeature);
        OGRFeature::DestroyFeature(pPondFeature);
        OGRFeature::DestroyFeature(pBratFeature);
    }

    pDem = (GDALDataset*) GDALOpen(inDem, GA_ReadOnly);
    int nRows, nCols;
    nCols = pDem->GetRasterXSize();
    nRows = pDem->GetRasterYSize();
    double transform[6];
    pDem->GetGeoTransform(transform);

    const char *waterRas = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/ponds_lidar.tif";
    const char *waterPts = "E:/etal/Projects/NonLoc/Beaver_Modeling/z_crap/pond_pts.shp";

    loadRasterValuesFromPoint(waterRas, waterPts, nRows, nCols, transform);

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

int createDamPoints(const char *bratPath, const char *pointPath)
{

}

int getRasterCol(double transform[6], double xCoord)
{
    int col;
    double xOffset, xDiv;

    xOffset = xCoord - transform[0];
    xDiv = xOffset/transform[1];
    col = floor(xDiv);

    return col;
}

int getRasterRow(double transform[6], double yCoord)
{
    int row;

    double yOffset, yDiv;

    yOffset = transform[3] - yCoord;
    yDiv = yOffset/fabs(transform[5]);
    row = floor(yDiv);

    return row;
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

    QString damElevName = "dam_elev";
    int nDamEl, nDemEl, row, col;
    QString demElevName = "GRID_CODE";
    QString fieldName;
    double damElev, demElev, watDep;

    OGRFeatureDefn *featDfn = layer->GetLayerDefn();

    for (int i=0; i<featDfn->GetFieldCount(); i++)
    {
        fieldName = QString::fromUtf8(featDfn->GetFieldDefn(i)->GetNameRef());
        if (fieldName == damElevName)
        {
            nDamEl = i;
        }
        else if (fieldName == demElevName)
        {
            nDemEl = i;
        }
    }

    OGRFeature *feature;

    float *depth = (float*) CPLMalloc(sizeof(float)*1);

    layer->ResetReading();
    while((feature = layer->GetNextFeature()) != NULL)
    {
        damElev = feature->GetFieldAsDouble(nDamEl);
        demElev = feature->GetFieldAsDouble(nDemEl);
        watDep = damElev - demElev;

        OGRPoint *point = (OGRPoint*) feature->GetGeometryRef();

        if (watDep > 0.0)
        {
            col = getRasterCol(transform, point->getX());
            row = getRasterRow(transform, point->getY());
            *depth = watDep;

            pRaster->GetRasterBand(1)->RasterIO(GF_Write, col, row, 1, 1, depth, 1, 1, GDT_Float32, 0,0);
        }
    }

    OGRFeature::DestroyFeature(feature);
    OGRDataSource::DestroyDataSource(pPoints);

    CPLFree(depth);

    GDALClose(pRaster);
}
