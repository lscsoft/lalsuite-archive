/*
 * $Id$
 *
 * Copyright (C) 2006  Kipp C. Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                         tokenizer.Tokenizer Class
 *
 * ============================================================================
 */


#include <Python.h>
#include <ctype.h>
#include <stdlib.h>
#include <tokenizer.h>


/*
 * ============================================================================
 *
 *                               Tokenizer Type
 *
 * ============================================================================
 */


#define STREAM_QUOTE '"'
#define STREAM_ESCAPE '\\'


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	/* list of the types to which parsed tokens will be converted */
	PyObject **types;
	/* end of the types list */
	PyObject **types_length;
	/* the type to which the next parsed token will be converted */
	PyObject **type;
	/* delimiter character to be used in parsing */
	Py_UNICODE delimiter;
	/* size of internal buffer, minus null terminator */
	int allocation;
	/* internal buffer */
	Py_UNICODE *data;
	/* end of internal buffer's contents (null terminator) */
	Py_UNICODE *length;
	/* current offset in buffer */
	Py_UNICODE *pos;
} ligolw_Tokenizer;


/*
 * Append the contents of a unicode object to a tokenizer's internal
 * buffer, increasing the size of the buffer if needed.
 */


static int add_to_data(ligolw_Tokenizer *tokenizer, PyObject *unicode)
{
	int n = PyUnicode_GET_SIZE(unicode);

	if(n) {
		if(tokenizer->length - tokenizer->data + n > tokenizer->allocation) {
			/*
			 * convert pointers to integer offsets
			 */

			int pos = tokenizer->pos - tokenizer->data;
			int length = tokenizer->length - tokenizer->data;

			/*
			 * increase buffer size, adding 1 to leave room for
			 * the null terminator
			 */

			Py_UNICODE *old_data = tokenizer->data;

			tokenizer->data = realloc(tokenizer->data, (tokenizer->allocation + n + 1) * sizeof(*tokenizer->data));
			if(!tokenizer->data) {
				/*
				 * memory failure, restore pointer and exit
				 */

				tokenizer->data = old_data;
				return -1;
			}
			tokenizer->allocation += n;

			/*
			 * convert integer offsets back to pointers
			 */

			tokenizer->pos = &tokenizer->data[pos];
			tokenizer->length = &tokenizer->data[length];
		}

		/*
		 * copy data from unicode into buffer, appending null
		 * terminator
		 */

		memcpy(tokenizer->length, PyUnicode_AsUnicode(unicode), n * sizeof(*tokenizer->length));
		tokenizer->length += n;
		*tokenizer->length = 0;
	}

	/*
	 * success
	 */

	return 0;
}


/*
 * Shift the contents of the tokenizer's buffer so that the data starting
 * at pos is moved to the start of the buffer.  When moving data, add 1 to
 * the length to also move the null terminator.
 */


static void advance_to_pos(ligolw_Tokenizer *tokenizer)
{
	if(tokenizer->pos != tokenizer->data) {
		tokenizer->length -= tokenizer->pos - tokenizer->data;
		memmove(tokenizer->data, tokenizer->pos, (tokenizer->length - tokenizer->data + 1) * sizeof(*tokenizer->data));
		tokenizer->pos = tokenizer->data;
	}
}


/*
 * Free the tokenizer's types list.
 */


static void unref_types(ligolw_Tokenizer *tokenizer)
{
	for(tokenizer->type = tokenizer->types; tokenizer->type < tokenizer->types_length; tokenizer->type++)
		Py_DECREF(*tokenizer->type);

	free(tokenizer->types);
	tokenizer->types = NULL;
	tokenizer->types_length = NULL;
	tokenizer->type = NULL;
}


/*
 * Identify the next token to extract from the tokenizer's internal buffer.
 * On success, start will be left pointing to the address of the start of
 * the string, and end will be pointing to the first character after the
 * string.  If no token is encountered, only whitespace between two
 * delimiters, then start and end are both set to NULL.  The return value
 * is the Python type to which the text should be converted, or NULL on
 * error.  On error, the values of start and end are undefined.  Raises
 * StopIteration if the end of the tokenizer's internal buffer is reached,
 * or ValueError if a parse error occurs.
 */


static PyObject *next_token(ligolw_Tokenizer *tokenizer, Py_UNICODE **start, Py_UNICODE **end)
{
	Py_UNICODE *pos = tokenizer->pos;
	Py_UNICODE *bailout = tokenizer->length;
	PyObject *type = *tokenizer->type;

	/*
	 * The following code matches the pattern:
	 *
	 * any amount of white-space + " + non-quote characters + " + any
	 * amount of white-space + delimiter
	 *
	 * or
	 *
	 * any amount of white-space + non-quote, non-white-space,
	 * non-delimiter characters + any amount of white-space + delimiter
	 *
	 * The middle bit is returned as the token.  '"' characters can be
	 * escaped by preceding them with a '\' character.
	 */

	/*
	 * start == a white-space to non-white-space transition outside of
	 * a quote, or a non-quoted to quoted transition.
	 *
	 * end == a non-white-space to white-space transition outside of a
	 * quote, or a delimiter outside of a quote, or a quoted to
	 * non-quoted transition.
	 */

	if(pos >= bailout)
		goto stop_iteration;
	while(Py_UNICODE_ISSPACE(*pos))
		if(++pos >= bailout)
			goto stop_iteration;
	if(*pos == STREAM_QUOTE) {
		int escaped = 0;
		*start = ++pos;
		if(pos >= bailout)
			goto stop_iteration;
		while((*pos != STREAM_QUOTE) || escaped) {
			escaped = (*pos == STREAM_ESCAPE) && !escaped;
			if(++pos >= bailout)
				goto stop_iteration;
		}
		*end = pos;
		if(++pos >= bailout)
			goto stop_iteration;
	} else {
		*start = pos;
		/* FIXME:  why the != STREAM_QUOTE? */
		while(!Py_UNICODE_ISSPACE(*pos) && (*pos != tokenizer->delimiter) && (*pos != STREAM_QUOTE))
			if(++pos >= bailout)
				goto stop_iteration;
		*end = pos;
		if(*start == *end)
			/* nothing but unquoted whitespace between
			 * delimiters */
			*start = *end = NULL;
	}
	while(*pos != tokenizer->delimiter) {
		if(!Py_UNICODE_ISSPACE(*pos))
			goto parse_error;
		if(++pos >= bailout)
			goto stop_iteration;
	}

	/*
	 * After this, tokenizer->pos points to the first character after
	 * the delimiter that terminated this current token.
	 */

	tokenizer->pos = ++pos;

	/*
	 * Select the next type
	 */

	if(++tokenizer->type >= tokenizer->types_length)
		tokenizer->type = tokenizer->types;

	/*
	 * Done.  *start points to the first character of the token, *end
	 * points to the first character following the token (or both are
	 * NULL if there was nothing but unquoted whitespace),
	 * tokenizer->pos and tokenizer->type have been advanced in
	 * readiness for the next token, and the return value is the python
	 * type to which the current token is to be converted.
	 */

	return type;

	/*
	 * Errors
	 */

stop_iteration:
	advance_to_pos(tokenizer);
	PyErr_SetNone(PyExc_StopIteration);
	return NULL;

parse_error:
	{
	PyObject *unicode = PyUnicode_FromUnicode(*start, tokenizer->length - *start);
	PyErr_SetObject(PyExc_ValueError, unicode);
	Py_DECREF(unicode);
	return NULL;
	}
}


/*
 * append() method
 */


static PyObject *append(PyObject *self, PyObject *data)
{
	int fail;

	if(PyUnicode_Check(data)) {
		fail = add_to_data((ligolw_Tokenizer *) self, data);
	} else if(PyString_Check(data)) {
		/* decode to unicode */
		if(!(data = PyUnicode_FromObject(data)))
			/* decode failure */
			return NULL;
		fail = add_to_data((ligolw_Tokenizer *) self, data);
		Py_DECREF(data);
	} else {
		PyErr_SetObject(PyExc_TypeError, data);
		return NULL;
	}

	if(fail < 0)
		return PyErr_NoMemory();

	Py_INCREF(self);
	return self;
}


/*
 * __del__() method
 */


static void __del__(PyObject *self)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;

	unref_types(tokenizer);
	free(tokenizer->data);
	tokenizer->data = NULL;
	tokenizer->allocation = 0;
	tokenizer->length = NULL;
	tokenizer->pos = NULL;

	self->ob_type->tp_free(self);
}


/*
 * __init__() method
 */


static int __init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	PyObject *arg;

	if(!PyArg_ParseTuple(args, "U", &arg))
		return -1;

	if(PyUnicode_GET_SIZE(arg) != 1) {
		PyErr_SetString(PyExc_ValueError, "len(delimiter) != 1");
		return -1;
	}

	tokenizer->delimiter = *PyUnicode_AS_UNICODE(arg);
	tokenizer->types = malloc(1 * sizeof(*tokenizer->types));
	tokenizer->types_length = &tokenizer->types[1];
	tokenizer->types[0] = (PyObject *) &PyUnicode_Type;
	Py_INCREF(tokenizer->types[0]);
	tokenizer->type = tokenizer->types;
	tokenizer->allocation = 0;
	tokenizer->data = NULL;
	tokenizer->length = tokenizer->data;
	tokenizer->pos = tokenizer->data;

	return 0;
}


/*
 * __iter__() method
 */


static PyObject *__iter__(PyObject *self)
{
	Py_INCREF(self);
	return self;
}


/*
 * next() method
 */


static PyObject *next(PyObject *self)
{
	PyObject *type;
	PyObject *token;
	Py_UNICODE *start, *end;

	/*
	 * Identify the start and end of the next token.
	 */

	do {
		type = next_token((ligolw_Tokenizer *) self, &start, &end);
		if(!type)
			return NULL;
	} while(type == Py_None);

	/*
	 * Null-terminate the token.
	 */

	if(end)
		*end = 0;

	/*
	 * Extract token as desired type.
	 */

	if(start == NULL) {
		/* unquoted zero-length string == None */
		Py_INCREF(Py_None);
		token = Py_None;
	} else if(type == (PyObject *) &PyFloat_Type) {
		char ascii_buffer[end - start + 1];
		char *ascii_end;
		if(PyUnicode_EncodeDecimal(start, end - start, ascii_buffer, NULL))
			return NULL;
		token = PyFloat_FromDouble(strtod(ascii_buffer, &ascii_end));
		if(ascii_end == ascii_buffer || *ascii_end != 0) {
			/* strtod() couldn't convert the token, emulate
			 * float()'s error message */
			Py_DECREF(token);
			PyErr_Format(PyExc_ValueError, "invalid literal for float(): '%s'", ascii_buffer);
			token = NULL;
		}
	} else if(type == (PyObject *) &PyUnicode_Type) {
		token = PyUnicode_FromUnicode(start, end - start);
	} else if(type == (PyObject *) &PyString_Type) {
		token = PyUnicode_Encode(start, end - start, NULL, NULL);
	} else if(type == (PyObject *) &PyInt_Type) {
		token = PyInt_FromUnicode(start, end - start, 0);
	} else if(type == (PyObject *) &PyLong_Type) {
		token = PyLong_FromUnicode(start, end - start, 0);
	} else {
		token = PyObject_CallFunction(type, "u#", start, end - start);
	}

	/*
	 * Done.
	 */

	return token;
}


/*
 * set_types() method
 */


static PyObject *set_types(PyObject *self, PyObject *sequence)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	int length, i;

	/*
	 * Simplify the sequence access.
	 */

	sequence = PySequence_Tuple(sequence);
	if(!sequence)
		return NULL;
	length = PyTuple_GET_SIZE(sequence);

	/*
	 * Free the current internal type list.
	 */

	unref_types(tokenizer);

	/*
	 * Allocate a new internal type list.
	 */

	tokenizer->types = malloc(length * sizeof(*tokenizer->types));
	if(!tokenizer->types) {
		Py_DECREF(sequence);
		return PyErr_NoMemory();
	}
	tokenizer->type = tokenizer->types;
	tokenizer->types_length = &tokenizer->types[length];

	/*
	 * Copy the input sequence's contents into the internal type list.
	 */

	for(i = 0; i < length; i++) {
		tokenizer->types[i] = PyTuple_GET_ITEM(sequence, i);
		Py_INCREF(tokenizer->types[i]);
	}

	/*
	 * Done.
	 */

	Py_DECREF(sequence);
	Py_INCREF(Py_None);
	return Py_None;
}


/*
 * Attribute access.
 */


static PyObject *attribute_get_data(PyObject *obj, void *data)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) obj;

	return PyUnicode_FromUnicode(tokenizer->data, tokenizer->length - tokenizer->data);
}


/*
 * Type information
 */


static struct PyMethodDef methods[] = {
	{"append", append, METH_O, "Append a unicode object to the tokenizer's internal buffer.  Also accepts str objects as input."},
	{"set_types", set_types, METH_O, "Set the types to be used cyclically for token parsing.  This function accepts an iterable of callables.  Each callable will be passed the token to be converted as a unicode string.  Special fast-paths are included to handle the Python builtin types float, int, long, str, and unicode.  The default is to return all tokens as unicode objects."},
	{NULL,}
};


static struct PyGetSetDef getset[] = {
	{"data", attribute_get_data, NULL, "The current contents of the internal buffer.", NULL},
	{NULL,}
};


PyTypeObject ligolw_Tokenizer_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_Tokenizer),
	.tp_dealloc = __del__,
	.tp_doc =
"A tokenizer for LIGO Light Weight XML Stream and Array elements.  Converts\n" \
"(usually comma-) delimited text streams into sequences of Python objects.  An\n" \
"instance is created by calling the class with the delimiter character as the\n" \
"single argument.  Text is appended to the internal buffer by passing it to the\n" \
"append() method.  Tokens are extracted by iterating over the instance.  The\n" \
"Tokenizer is able to directly extract tokens as various Python types.  The\n" \
"set_types() method is passed a sequence of the types to which tokens are to be\n" \
"converted.  The types will be used in order, cyclically.  For example, passing\n" \
"[int] to set_types() causes all tokens to be converted to ints, while\n" \
"[str, int] causes the first token to be returned as a string, the second as an\n" \
"int, then the third as a string again, and so on.  The default is to extract\n" \
"all tokens as strings.\n" \
"\n" \
"Example:\n" \
"\n" \
">>> from glue.ligolw import tokenizer\n" \
">>> t = tokenizer.Tokenizer(\",\")\n" \
">>> t.set_types([str, int])\n" \
">>> list(t.append(\"a,10,b,2\"))\n" \
"['a', 10, 'b']\n" \
">>> list(t.append(\"0,\"))\n" \
"[20]\n" \
"\n" \
"Notes.  The last token will not be extracted until a delimiter character is\n" \
"seen to terminate it.  Tokens can be quoted with '\"' characters, which will\n" \
"removed before conversion to the target type.  An empty token (two delimiters\n" \
"with only whitespace between them) is returned as None regardless of the\n" \
"requested type.  To prevent a zero-length string token from being interpreted\n" \
"as None, place it in quotes.",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
	.tp_init = __init__,
	.tp_iter = __iter__,
	.tp_iternext = next,
	.tp_getset = getset,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".Tokenizer",
	.tp_new = PyType_GenericNew,
};
