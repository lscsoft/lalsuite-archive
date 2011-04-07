//lalinference.h

int setLALIFODataFromData(LALIFOData* internal,PyObject* old,PyObject* new);

PyObject* getLALIFODataFromData(LALIFOData* internal,PyObject* owner);
int setBarycenterInputFromData(BarycenterInput* internal,PyObject* old,PyObject* new);

PyObject* getBarycenterInputFromData(BarycenterInput* internal,PyObject* owner);

int setREAL8WindowFromData(REAL8Window* internal,PyObject* old,PyObject* new);

PyObject* getREAL8WindowFromData(REAL8Window* internal,PyObject* owner);

int setREAL8FFTPlanFromData(REAL8FFTPlan* internal,PyObject* old,PyObject* new);

PyObject* getREAL8FFTPlanFromData(REAL8FFTPlan* internal,PyObject* owner);

int setCOMPLEX16TimeSeriesFromData(COMPLEX16TimeSeries* internal,PyObject* old,PyObject* new);

PyObject* getCOMPLEX16TimeSeriesFromData(COMPLEX16TimeSeries* internal,PyObject* owner);

int setCOMPLEX16FrequencySeriesFromData(COMPLEX16FrequencySeries* internal,PyObject* old,PyObject* new);

PyObject* getCOMPLEX16FrequencySeriesFromData(COMPLEX16FrequencySeries* internal,PyObject* owner);

int setREAL8FrequencySeriesFromData(REAL8FrequencySeries* internal,PyObject* old,PyObject* new);

PyObject* getREAL8FrequencySeriesFromData(REAL8FrequencySeries* internal,PyObject* owner);

int setREAL8TimeSeriesFromData(REAL8TimeSeries* internal,PyObject* old,PyObject* new);

PyObject* getREAL8TimeSeriesFromData(REAL8TimeSeries* internal,PyObject* owner);

int setLALDetectorFromData(LALDetector* target,PyObject* old,PyObject* origin);

PyObject* getLALDetectorFromData(LALDetector internal);

int setLIGOTimeGPSFromData(LIGOTimeGPS target,PyObject* origin);
