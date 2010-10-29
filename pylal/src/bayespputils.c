/*
 * ============================================================================
 *
 *                C extensions for Bayesian post-processing codes
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>

#define MODULE_NAME "_bayespputils"

/* doc string */
const char BUDocstring[] =
"This module provides C extensions for Bayesian analysis and post-\n"
"processing codes.";

/* calculateConfidenceLevels
 * C extension function replacing pylal.bayespputils.calculateConfidenceLevels_slow .
 */
static PyObject *_calculate_confidence_levels(PyObject *self, PyObject *args)
{
    /***DECLARATION***/

    //Input objects
    PyArrayObject *shist_Array=NULL;
    PyArrayObject *skypoints_Array=NULL;
    PyObject* injbin;//bin in which the injection was found (if any)
    double skyres[1];
    int Nsamples[1];

    //internal variables
    double frac=0.; //fraction of points binned
    int Nbins=0;//total number of bins used from skybins
    int maxpos=0;
    int maxbin=0; //stores maxbin through loop
    double injectionconfidence=-0.1;
    int i; //iterators
    PyObject* py_new_toppoints_tuple=NULL;

    //Return values
    PyObject* py_injectionconfidence=NULL;
    PyObject* toppoints_list=NULL;

    //Create python lists:
    toppoints_list=PyList_New((Py_ssize_t)0);

    /***PARSE/PROCESS INPUT***/

    if (!PyArg_ParseTuple(args,"O!O!Odi",&PyArray_Type, &shist_Array,&PyArray_Type, &skypoints_Array,&injbin,skyres,Nsamples))  return NULL;

    //array iterators
    int shist_Array_stride=shist_Array->strides[0];

    int lenbins = shist_Array->dimensions[0];


    int skypoints_Array_stride_i=skypoints_Array->strides[0];
    int skypoints_Array_stride_j=skypoints_Array->strides[1];

    int shist_Array_current_position;
    int shist_Array_current_value;

    int skypoints_Array_current_position;

    int maxbin_count=0;

	/***MAIN LOOP
	 * This will look over the sample bins repeatedly determining the bin containing
	 * the highest number of samples and adding this to the fractional total up to
	 * the current confidence level.
	 * ***/
	while(frac<0.9999){
		maxbin=0;
		//Find bin with highest count...
		for(i=0;i<lenbins;i++){
			shist_Array_current_position=i*shist_Array_stride;
			shist_Array_current_value=*(int*)(shist_Array->data+shist_Array_current_position);
			if(shist_Array_current_value>maxbin){
				maxbin=shist_Array_current_value;
				maxpos=i;
			}
		}
					//...add this to the total and zero the bin.
		shist_Array_current_position=maxpos*shist_Array_stride;
		*(int*)(shist_Array->data+shist_Array_current_position)=0;
		maxbin_count+=maxbin;
		frac = frac + ( ((double)maxbin) / ((double)Nsamples[0]) ) ;
		Nbins++;

		//Build new toppoint tuple
		skypoints_Array_current_position=skypoints_Array_stride_i*maxpos;

		py_new_toppoints_tuple=Py_BuildValue("(ddid)",*(double*)(skypoints_Array->data+skypoints_Array_current_position),*(double*)(skypoints_Array->data+skypoints_Array_current_position+skypoints_Array_stride_j),maxpos,frac);

		//...add to toppoints list
		PyList_Append(toppoints_list,py_new_toppoints_tuple);
		Py_DECREF(py_new_toppoints_tuple);

		//If an injection bin was given see if this is in the current maxbin:
		Py_INCREF(Py_None);
		if(injbin!=Py_None){
			long int injbin_l = PyInt_AsLong(injbin);
			if(injbin_l == (long int)maxpos){
				injectionconfidence=frac;
				printf("Injection point found at confidence level %lf.\n",injectionconfidence);
			}
		}
		Py_DECREF(Py_None);
	}
	/**END MAIN LOOP**/

    //Create PyFloat/PyNone of injection confidence variable

    if( injectionconfidence > 0.) py_injectionconfidence=PyFloat_FromDouble(injectionconfidence);
    else {
        Py_INCREF(Py_None);
        py_injectionconfidence=Py_None;
    }

    //return value
    return Py_BuildValue("(OO)",py_injectionconfidence,toppoints_list);

}

/* skyhist_cart
 * C extension function replacing pylal.bayespputils.skyhist_cart .
 */
static PyObject *_skyhist_cart(PyObject *self, PyObject *args)
{
    /***DECLARATION***/

    //INPUT
    PyArrayObject *SkyCartsArray=NULL,*SamplesArray=NULL;

    //OUTPUT
    PyArrayObject *BinsArray=NULL;

    //INTERNAL VARIABLES
    int skybin_idx,sample_no;//iterators
    int a,b;//SkyCarts dimensions
    int m,n;//Samples dimensions

    //containers for loop variables
    double xsample,ysample,zsample,xskybin,yskybin,zskybin;
    double longi,lat,maxbin_value,dot;
    int maxbin;

    /***PARSE/PROCESS INPUT***/

    if (!PyArg_ParseTuple(args, "O!O!",&PyArray_Type, &SkyCartsArray,&PyArray_Type, &SamplesArray))  return NULL;

    if(SkyCartsArray->nd != 2 || SkyCartsArray->descr->type_num != PyArray_DOUBLE) {
        PyErr_SetObject(PyExc_ValueError, (PyObject *) SkyCartsArray);
        return NULL;
    }

    if( (SamplesArray->nd !=2 ) || SamplesArray->descr->type_num != PyArray_DOUBLE) {
        PyErr_SetObject(PyExc_ValueError, (PyObject *) SamplesArray);
        return NULL;
    }

    //get dimensions
    a=SkyCartsArray->dimensions[0];
    b=SkyCartsArray->dimensions[1];

    m=SamplesArray->dimensions[0];
    n=SamplesArray->dimensions[1];


    npy_intp bin_dims[1] = {a};

    //create bins vector
    BinsArray=(PyArrayObject *) PyArray_SimpleNew(1,bin_dims,PyArray_INT);
    Py_INCREF(BinsArray);

    if(BinsArray==NULL) return NULL;

    //get iterators
    int SamplesArray_stride_i=SamplesArray->strides[0];
    int SamplesArray_stride_j=SamplesArray->strides[1];

    int SkyCartsArray_stride_i=SkyCartsArray->strides[0];
    int SkyCartsArray_stride_j=SkyCartsArray->strides[1];

    int BinsArray_stride_i=BinsArray->strides[0];

    int SkyCartsArray_iterx;

    //zero the bins
    for(skybin_idx=0;skybin_idx<a;skybin_idx++){
        *(int*)(BinsArray->data+skybin_idx*BinsArray_stride_i) = 0;
    }

    int SamplesArray_stride_j_RA_iter=0;
    int SamplesArray_stride_j_dec_iter=SamplesArray_stride_j;


    /***MAIN LOOP***/
    /* Now going to loop over samples. For each sapmle the max. dot product
     * corresponds to the sky bin that its going to be associated with.
     */

    for(sample_no=0;sample_no<m;sample_no++){

        longi = *(double*)(SamplesArray->data+sample_no*SamplesArray_stride_i+SamplesArray_stride_j_RA_iter);
        lat = *(double*)(SamplesArray->data+sample_no*SamplesArray_stride_i+SamplesArray_stride_j_dec_iter);

        //Cartesian co of sample for dot product
        xsample=cos(lat)*cos(longi);
        ysample=cos(lat)*sin(longi);
        zsample=sin(lat);
        maxbin=0;
        maxbin_value=-1.;
        for(skybin_idx=0;skybin_idx<a;skybin_idx++){

            SkyCartsArray_iterx=skybin_idx*SkyCartsArray_stride_i;

            //Cartesian co of bin for dot product
            xskybin = *(double*)(SkyCartsArray->data+SkyCartsArray_iterx);
            yskybin = *(double*)(SkyCartsArray->data+SkyCartsArray_iterx+SkyCartsArray_stride_j);
            zskybin = *(double*)(SkyCartsArray->data+SkyCartsArray_iterx+2*SkyCartsArray_stride_j);

            //dot (x,y,z) with each skybin
            dot=( xsample*xskybin + ysample*yskybin + zsample*zskybin );
            if(  dot > maxbin_value ){//if dot is the biggest so far update the maxbin info
                maxbin = skybin_idx;
                maxbin_value = dot;

            }

        }

        //add one to the max bin count for the highest scoring dot product
        *(int*)(BinsArray->data+maxbin*BinsArray_stride_i) += 1;
    }

    return PyArray_Return(BinsArray);

}

/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef BUmethods[] = {
    {"_skyhist_cart", _skyhist_cart, METH_VARARGS,
    "This function bins samples on a grid generated from pylal.skylocutils\n"
    "from the marginal posterior for the sky position by finding the maximum\n"
    "dot product from each sky bin for each sample."
    },
    {"_calculate_confidence_levels", _calculate_confidence_levels, METH_VARARGS,
    "Given a histogram of posterior samples this will determine confidence levels"
    "using a greedy algorithm."
    },
    {NULL,}
};


void init_bayespputils(void)
{
    (void)Py_InitModule3(MODULE_NAME,BUmethods,BUDocstring);
    import_array();
}
