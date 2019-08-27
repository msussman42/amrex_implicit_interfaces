//  December 30, 2003
//  TO USE IN ANOTHER CODE: Check carefully use of copyTo and createSubArray()
//
//----------------------------------------------------------------------
// May 2, 2002
// Change definition of ArrayT (with templates) so that constructor
//    is ArrayT(npts, origin) where npts is tot nb of points in each dimension
//    For compatibility with current version, create new constructor
//    ArrayT(Vec3i& npts, Vec3i origin, T* array=0);

// When going over loop indices: for (i=getOrigin()[0]; i < getMaxDims()[0]; i++) {, ...

#ifndef _ARRAYT_H_
#define _ARRAYT_H_

#include <stdio.h>
#include "Vec3i.h"

template <class T>
class ArrayT {
private:
    T* data; // data allocated in constructor, unless float* pointer passed to constructor
    int npts; // total number of points in the array
    int dims[3];
    int n1, n2, n3;  // first index varies fastest (like Fortran)
    int np1, np2, np3;  // tot nb points in each dimension
    int n1m, n2m, n3m;  // 
    int externalPtr;
    Vec3i maxDims;
    Vec3i origin;
    Vec3i totPts; // maxDims - origin
public:
    enum {CSTYLE=0,FORTRANSTYLE};
    // allocation, Fortran style
    // args: max1, max2, max3, min1=min2=min3=0
	ArrayT() {printf("ArrayT:: default constructor\n");}
    ArrayT(T* array, const Vec3i& n_);
    ArrayT(T* array, const Vec3i& n_, const Vec3i& nm_);
    ArrayT(const Vec3i& n_);
    ArrayT(const Vec3i& n_, const Vec3i& nm_);
    ArrayT(T* array, const int* n_) {
        init(array, n_[0], n_[1], n_[2]);
    }
    ArrayT(T* array, int n1_, int n2_, int n3_) {
    // Attach ArrayT object to external pointer
        init(array, n1_, n2_, n3_);
    }
    ArrayT(int n1_, int n2_, int n3_, int n1m_=0, int n2m_=0, int n3m_=0) {
        init(0, n1_, n2_, n3_, n1m_, n2m_, n3m_); 
    }
    ~ArrayT() {
        remove();
    }
    inline void remove() {
        if (externalPtr == 0) {
            delete [] data;
        }
    }
    // ATTENTION: init() changed to zero default
    void init(T* array, int n1_, int n2_, int n3_, int n1m_=0, int n2m_=0, int n3m_=0) {
        n1m = n1m_; // min dimension value, origin
        n2m = n2m_;
        n3m = n3m_;

        origin.setValue(n1m, n2m, n3m);

        np1 = n1_; // size in dimension 1
        np2 = n2_;
        np3 = n3_;

        totPts.setValue(np1, np2, np3);

        // maximum dimension: array[n1m:n1, n2m:n2, n3m:n3]
        n1 = n1m + np1; // C-style: origin[0], origin[0]+1, ...,np1-1
        n2 = n2m + np2; // C-style: origin[0], origin[0]+1, ...,np1-1
        n3 = n3m + np3; // C-style: origin[0], origin[0]+1, ...,np1-1

        maxDims = origin + totPts;
        npts = np1*np2*np3;

        if (array) {
            externalPtr = 1;
            data = array;
        } else {
            externalPtr = 0;
            data = new T [npts];
        }
    }
    int* getDims() {dims[0] = np1; dims[1] = np2; dims[2] = np3; return dims;}
    Vec3i& getMaxDims() {return maxDims;}
    inline T* getDataPtr() {return data;}
    inline int getSize() {return npts;}
    inline int getSize() const {return npts;}
	// NOT SURE: T or T& ???
    inline T& get(int i, int j=0, int k=0) {
        return data[(i-n1m)+np1*((j-n2m)+np2*(k-n3m))];
    }
	// const arguments required
    inline T& get(const int i, const int j=0, const int k=0) const {
		printf("const get\n");
        return data[(i-n1m)+np1*((j-n2m)+np2*(k-n3m))];
    }
    //inline void set(int i, int j, int k, const T& value) {
    void set(int i, int j, int k, const T& value) {
        data[(i-n1m)+np1*((j-n2m)+np2*(k-n3m))] = value;
    }
	ArrayT<T>& operator=(const  ArrayT<T>& arr) const;
	ArrayT<T>& operator=(ArrayT<T>& arr);
    inline T operator[](int i) { return data[i]; }
	T& operator()(int i, int j=0, int k=0) { return get(i,j,k);}
	T& operator()(int i, int j=0, int k=0) const { return get(i,j,k);}
// Copy from this to array
    int copyTo(ArrayT& array);
    Vec3i& getOrigin() {return origin;} // expensive (2 objects)
    Vec3i& getTotPoints() {return totPts;}
// Copy to this from array
    int copyFrom(ArrayT& array);
// Copy from this to array
    //int copyTo(ArrayT& array, Vec3i& toOrig, Vec3i fromOrig, Vec3i fromRange);
    int copyTo(ArrayT& array, const Vec3i& toOrig, Vec3i fromRange, Vec3i fromOrigin);
// Copy from array array
    int copyFrom(ArrayT& array, Vec3i fromOrig, Vec3i fromRange, const Vec3i& toOrigin);
// set to value
    void setTo(const T& value);
    void setOrigin(const Vec3i& origin) {
        n1m = origin[0];
        n2m = origin[1];
        n3m = origin[2];
        this->origin.setValue(n1m, n2m, n3m);
    }
	void printInfo(const char* msg) {
		//printf("--- %s ---\n", msg);
		getTotPoints().print(msg);
	}
    void print(char*);
    void printi(char*);
    void print(char*, Vec3i orig, Vec3i size);
    void printi(char*, Vec3i orig, Vec3i size);
    // can only resize arrays allocated internally
    void resize(Vec3i& arraySize, Vec3i& origin); // max nb points in 3 directions, origin
    void resize(Vec3i& arraySize); // max nb points in 3 directions, origin
    int getCapacity() {return npts;} // return total size of array (product of 3 dimensions) 

    // User is responsible for deleting allocated memory of subArray
	// copy by default
    ArrayT* createSubArray(Vec3i subOrigin, Vec3i subWidth, Vec3i newOrigin, int copy=1);
	ArrayT* createArray(); // no initialization
	ArrayT* createCopy();

	// GE: Dec. 10, 2003
	T maximum(); 
	T minimum();
	T minval(const Vec3i& min, const Vec3i& max);
	T maxval(const Vec3i& min, const Vec3i& max);


    //ArrayT& operator+(const ArrayT& a);
};


//----------------------------------------------------------------------
template <class T>
ArrayT<T>& ArrayT<T>::operator=(const ArrayT& arr) const {
	Vec3i t = arr.getTotPoints();
	Vec3i o = arr.getOrigin();
	// Always make a copy of the array (even if data in this points
	// to external data
	
	if (&arr != this) {
		init(0, t[0], t[1], t[2], o[0], o[1], o[2]);
		arr.copyTo(*this);
	}
	return *this;
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>& ArrayT<T>::operator=(ArrayT& arr) {
	Vec3i t = arr.getTotPoints();
	Vec3i o = arr.getOrigin();
	// Always make a copy of the array (even if data in this points
	// to external data
	
	if (&arr != this) {
		init(0, t[0], t[1], t[2], o[0], o[1], o[2]);
		arr.copyTo(*this);
	}
	return *this;
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>::ArrayT(T* array, const Vec3i& n_)
{
    init(array, n_[0], n_[1], n_[2]);
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>::ArrayT(T* array, const Vec3i& n_, const Vec3i& nm_)
{
    init(array, n_[0], n_[1], n_[2], nm_[0], nm_[1], nm_[2]);
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>::ArrayT(const Vec3i& n_)
{
    init(0, n_[0], n_[1], n_[2]);
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>::ArrayT(const Vec3i& n_, const Vec3i& nm_)
{
    init(0, n_[0], n_[1], n_[2], nm_[0], nm_[1], nm_[2]);
}
//----------------------------------------------------------------------
#if 0
template <class T>
void ArrayT<T>::computeMinMax()
{
    min = max = data[0];
    for (int i=0; i < npts; i++) {
        //printf("i,data= %d, %g, min/max= %g, %g\n", i, data[i], min, max);
        min = data[i] < min ? data[i] : min;
        max = data[i] > max ? data[i] : max;
    }
}
#endif
//----------------------------------------------------------------------
template <class T>
int ArrayT<T>::copyTo(ArrayT& a)
// copy to "a" from "this"
// return 0 if 0k, -1 otherwise
{
    int* dims = a.getDims();
    if (dims[0] != np1 || dims[1] != np2 || dims[2] != np3) {
        return(-1);
    }
	int npts = a.getSize();
    T* ptr = a.getDataPtr();
	// memcpy only works on primitive types. Need more general copy operation 
	// for general types
    //memcpy(ptr, data, sizeof(T)*npts);

	// more inefficient, but general. Requires that T have an assignment operator
	for (int i=0; i < npts; i++) {
		ptr[i] = data[i];
	}

    return(0);
}
//----------------------------------------------------------------------
template <class T>
int ArrayT<T>::copyFrom(ArrayT& a)
// copy from "a" to "this"
{
    //printf("enter copyFrom\n");
    int* dims = a.getDims();
    //printf("dims= %d, %d, %d\n", dims[0], dims[1], dims[2]);
    if (dims[0] != np1 || dims[1] != np2 || dims[2] != np3) {
        //theMsg->printf("error in copyFrom (array3D.cpp)\n");
        return(-1);
    }
    T* ptr = a.getDataPtr();
    memcpy(data, ptr, sizeof(T)*npts);
    return(0);
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::setTo(const T& value)
{
    for (int i=0; i < npts; i++) {
        data[i] = value;
    }
    //memset(data, (int) value, sizeof(GE_FLOAT)*npts);
}
//----------------------------------------------------------------------
template <class T>
int ArrayT<T>::copyTo(ArrayT& a, const Vec3i& toOrigin, Vec3i fromRange, Vec3i fromOrigin)
// Must pass last two arguments by value because they could otherwise change
//int ArrayT<T>::copyTo(ArrayT& a, const Vec3i& toOrigin, Vec3i& fromRange, Vec3i& fromOrigin)
// copy from "this" to "a"
// fromRange cannot be const: it is changed 4 lines down.
// fromOrigin cannot be const: a non-const function is applied to it
{
    Vec3i fromMaxDims = fromOrigin + fromRange;
    fromMaxDims.clampMaxTo(maxDims);
    fromOrigin.clampMinTo(origin);
    fromRange = fromMaxDims - fromOrigin;
    Vec3i toMaxDims = toOrigin + fromRange;

    if (toOrigin < a.getOrigin()) return -1;
    if ((toOrigin+fromRange) > a.getMaxDims()) return -1;

    for (int k=0; k < fromRange.z(); k++) {
    for (int j=0; j < fromRange.y(); j++) {
    for (int i=0; i < fromRange.x(); i++) {
        a.set(i+toOrigin.x(), j+toOrigin.y(), k+toOrigin.z(), 
                get(i+fromOrigin.x(), j+fromOrigin.y(), k+fromOrigin.z()));
    }}}
    return(0);
}
//----------------------------------------------------------------------
template <class T>
int ArrayT<T>::copyFrom(ArrayT& a, Vec3i fromOrig, Vec3i fromRange, const Vec3i& toOrigin)
{
    Vec3i fromMaxDims = fromOrig + fromRange;
    fromMaxDims.clampMaxTo(a.getMaxDims());
    fromOrig.clampMinTo(a.getOrigin());
    fromRange = fromMaxDims - fromOrig;
    Vec3i toMaxDims = toOrigin + fromRange;

    if (toOrigin < this->origin) return -1;
    if ((toOrigin+fromRange) > maxDims) return -1;

    for (int k=0; k < fromRange[2]; k++) {
    for (int j=0; j < fromRange[1]; j++) {
    for (int i=0; i < fromRange[0]; i++) {
        set(i+toOrigin.x(), j+toOrigin.y(), k+toOrigin.z(),  a.get(i+fromOrig.x(), j+fromOrig.y(), k+fromOrig.z()));
    }}}
    return(0);
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::resize(Vec3i& arraySize, Vec3i& origin) // max nb points in 3 directions, origin
{
    if (externalPtr) {
        printf("Array3d: Error: cannot resize external array\n");
        return;
    }
    if (data) {
        delete [] data;
        data = 0;
    }
    init(0, arraySize[0], arraySize[1], arraySize[2], origin[0], origin[1], origin[2]);
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::resize(Vec3i& arraySize) // max nb points in 3 directions, origin
{
    Vec3i zero = Vec3i(0,0,0);
    resize(arraySize, zero); // zero only gets deleted after return from resize, so ok
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>* ArrayT<T>::createArray()
{
    return( new ArrayT(getTotPoints(), getOrigin()) );
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>* ArrayT<T>::createCopy()
{
    ArrayT<T>* arr = new ArrayT(getTotPoints(), getOrigin());
	copyTo(*arr);
	return arr;
}
//----------------------------------------------------------------------
template <class T>
ArrayT<T>* ArrayT<T>::createSubArray(Vec3i subOrigin, Vec3i subWidth, Vec3i newOrigin, int copy)
{
// Given the current array, cut out an array of size subWidth at origin subOrigin.
// Create new array set to the cutout with an origin at newOrigin
// The method returns new origin and width of subarray in subOriginand subWidth
// Copy contents of original array into subarray (perhaps put an options to prevent copying?)
//
// Do not allow arguments to change. Why? Because sometimes the user makes use of
// values from the class object, and these should not change!

    // keep origin within bounds
    Vec3i orig = subOrigin;
    Vec3i maxDimens = subOrigin + subWidth;
	//maxDimens.print("maxDimens: ");
	//subOrigin.print("subOrigin: ");
    if (subOrigin > maxDims) return 0;

    orig.clampMinTo(origin);
	//totPts.print("totPts: ");
    maxDimens.clampMaxTo(maxDims);
	//maxDimens.print("maxDimens: ");

	//orig.print("orig: ");

    //Vec3i size = maxDimen - subOrigin;
    subOrigin = orig;
    subWidth = maxDimens - orig;

	//subWidth.print("newOrigin: ");
	//subWidth.print("subOrigin: ");
	//subWidth.print("subWidth: ");

	ArrayT<T>* arr = new ArrayT(subWidth, newOrigin);

	if (copy) {
		copyTo(*arr, newOrigin, subWidth, subOrigin);
	}

    return(arr);
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::printi(char* msg)
{
	printi(msg, getOrigin(), getTotPoints());
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::print(char* msg)
{
	print(msg, getOrigin(), getTotPoints());
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::printi(char* msg, Vec3i orig, Vec3i size)
{
	// size: number of elements to print

	//Vec3i mx = orig + size;
	Vec3i mx = orig + size;
	size.print("size: ");
	orig.print("orig: ");
	mx.print("mx: ");
    printf("\n-------------- %s ---------------\n", msg);

    for (int k=orig[2]; k < mx[2]; k++) {
    for (int j=orig[1]; j < mx[1]; j++) {
    for (int i=orig[0]; i < mx[0]; i++) {
        printf("  i,j,k=  %d, %d, %d, arr= %d\n", i,j,k,get(i,j,k));
    }}}
    printf("----------------------------------\n");
}
//----------------------------------------------------------------------
template <class T>
void ArrayT<T>::print(char* msg, Vec3i orig, Vec3i size)
{
	// size: number of elements to print

	//Vec3i mx = orig + size;
	Vec3i mx = orig + size;
	size.print("size: ");
	orig.print("orig: ");
	mx.print("mx: ");
    printf("\n-------------- %s ---------------\n", msg);

    for (int k=orig[2]; k < mx[2]; k++) {
    for (int j=orig[1]; j < mx[1]; j++) {
    for (int i=orig[0]; i < mx[0]; i++) {
        printf("  i,j,k=  %d, %d, %d, arr= %g\n", i,j,k,get(i,j,k));
    }}}
    printf("----------------------------------\n");
}
//----------------------------------------------------------------------
template <class T>
T ArrayT<T>::maximum()
{
	T mx = data[0];
	for (int i=1; i < npts; i++) {
		mx = data[i] > mx ? data[i] : mx;
	}
	return mx;
}
//----------------------------------------------------------------------
template <class T>
T ArrayT<T>::minimum()
{
	T mn = data[0];
	for (int i=1; i < npts; i++) {
		mn = data[i] < mn ? data[i] : mn;
	}
	return mn;
}
//----------------------------------------------------------------------
template <class T>
T ArrayT<T>::minval(const Vec3i& min, const Vec3i& max)
// bounds are included in the search
{
// could be made more efficient

	int i, j, k;
	T mn = get(min[0], min[1], min[2]);
	for (int k=min[0]; k <= max[0]; k++) {
	for (int j=min[1]; j <= max[1]; j++) {
	for (int i=min[2]; i <= max[2]; i++) {
		mn = mn < get(i,j,k) ? get(i,j,k) : mn;
	}}}
	return mn;
}
//----------------------------------------------------------------------
template <class T>
T ArrayT<T>::maxval(const Vec3i& min, const Vec3i& max)
// bounds are included in the search
{
// could be made more efficient

	int i, j, k;
	T mx = get(min[0], min[1], min[2]);
	for (k=min[2]; k <= max[2]; k++) {
	for (j=min[1]; j <= max[1]; j++) {
	for (i=min[0]; i <= max[0]; i++) {
		mx = mx < get(i,j,k) ? get(i,j,k) : mx;
	}}}
	return mx;
}
//----------------------------------------------------------------------
/**
template <class T>
ArrayT<T>& ArrayT<T>::operator+(const ArrayT<T>& a)   // a is not to be modified
{
	int nbPts = this->getSize();
	int n     = a.getSize();
	if (n != nbPts) {
		printf("inconsistent array sizes in Array3D::operator+\n");
	}
	Vec3i& v = this->getTotPoints();
	ArrayT<T>* c = new ArrayT<T>(a.getTotPoints()); // returns Vec3i
	const T* af = a.getDataPtr();
	const T* bf = this->getDataPtr(); // (or simply data)
	T* cf = c->getDataPtr();
	for (int i=0; i < nbPts; i++) {
		cf[i] = af[i] + bf[i];
	}
	return *c;
}
***/
//----------------------------------------------------------------------
#ifdef STANDALONE
void main()
{
    float* x = new float [1000];
    ArrayT arr(x, 10, 10, 10);
    printf("arr->getDataPtr() = %ld\n", (long) arr.getDataPtr());
    exit(0);
}
#endif

// #ifndef _ARRAYT_H_
#endif
