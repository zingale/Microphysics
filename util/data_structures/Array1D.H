template<typename DataType, long size_i>
class Array1D {
public:
  DataType data[size_i];

  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  DataType& operator()(const int i = 1)
  {
    return data[i-1];
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  DataType* dataPtr()
  {
    return &data[0];
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  constexpr long size()
  {
    return size_i;
  }
};
