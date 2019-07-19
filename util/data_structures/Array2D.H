template<typename DataType, long size_i, long size_j>
class Array2D {
public:
  DataType data[size_i * size_j];

  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  DataType& operator()(const int i = 1, const int j = 1)
  {
    return data[(i - 1) + (j - 1) * size_i];
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  DataType* dataPtr()
  {
    return &data[0];
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  constexpr long size()
  {
    return size_i * size_j;
  }
};
