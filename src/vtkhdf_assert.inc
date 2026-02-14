
#ifdef NDEBUG
# define ASSERT(x) !! assert( x )
#else
# define ASSERT(x) if(.not.(x)) call vtkhdf_assert(__FILE__,__LINE__)
#endif

#define INSIST(x) if(.not.(x)) call vtkhdf_assert(__FILE__,__LINE__)
