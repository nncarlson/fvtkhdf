/* H5 interface functions that must be done from the C side */

#include "hdf5.h"

/* Functions to return object IDs that are created at run time */
hid_t H5P_DATASET_CREATE_value() { return H5P_DATASET_CREATE; }
hid_t H5P_GROUP_CREATE_value() { return H5P_GROUP_CREATE; }
hid_t H5P_FILE_ACCESS_value() {return H5P_FILE_ACCESS; }
hid_t H5P_CRT_ORDER_TRACKED_value() { return H5P_CRT_ORDER_TRACKED; }
hid_t H5P_CRT_ORDER_INDEXED_value() { return H5P_CRT_ORDER_INDEXED; }
hid_t H5P_DATASET_XFER_value() { return H5P_DATASET_XFER; }
hid_t H5T_NATIVE_UINT8_value() { return H5Tcopy(H5T_NATIVE_UINT8); }
hid_t H5T_NATIVE_INT32_value() { return H5Tcopy(H5T_NATIVE_INT32); }
hid_t H5T_NATIVE_INT64_value() { return H5Tcopy(H5T_NATIVE_INT64); }
hid_t H5T_NATIVE_FLOAT_value() { return H5Tcopy(H5T_NATIVE_FLOAT); }
hid_t H5T_NATIVE_DOUBLE_value() { return H5Tcopy(H5T_NATIVE_DOUBLE); }
hid_t H5T_NATIVE_CHARACTER_value() {
  hid_t type_id = H5Tcopy(H5T_FORTRAN_S1);
  H5Tset_size(type_id, 1);
  H5Tset_strpad(type_id, H5T_STR_SPACEPAD);
  return type_id;
}

#ifdef USE_MPI
#include "mpi.h"

/* Wrapper function that translates the input Fortran comm handle
   to a C handle, and hardwires MPI_INFO_NULL for info */
herr_t
H5Pset_fapl_mpio_Fcomm(hid_t fapl_id, MPI_Fint Fcomm) {
  MPI_Comm comm = MPI_Comm_f2c(Fcomm);
  return H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
}

/* Wrapper function that translates the returned C comm handle to
   a Fortran handle and skips getting the info handle entirely */
herr_t
H5Pget_fapl_mpio_Fcomm(hid_t fapl, MPI_Fint *Fcomm) {
  MPI_Comm comm;
  herr_t herr = H5Pget_fapl_mpio(fapl, &comm, NULL);
  if (herr < 0) comm = MPI_COMM_NULL;
  *Fcomm = MPI_Comm_c2f(comm);
  return herr;
}
#endif
