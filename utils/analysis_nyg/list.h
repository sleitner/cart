#ifndef __NG_LIST_H__
#define __NG_LIST_H__


namespace ng
{
  template<class T> class List
  {
  public:

    inline const T* Data() const { return p.List; }
    inline int Size() const { return p.Size; }

    List()
    {
      p.Size = p.BufferSize = 0;
      p.List = 0;
    }

    virtual ~List()
    {
      if(p.List != 0) delete [] p.List;
    }

    void Add(const T& w)
    {
      if(this->IsValid(w))
	{
	  if(p.Size == p.BufferSize)
	    {
	      p.BufferSize += 100;
	      T *tmp = new T[p.BufferSize];
	      cart_assert(tmp != 0);

	      int i;
	      for(i=0; i<p.Size; i++)
		{
		  tmp[i] = p.List[i];
		}
	      delete [] p.List;
	      p.List = tmp;
	    }

	  p.List[p.Size++] = w;
	}
      }

    void operator+=(const T& w)
    {
      this->Add(w);
    }

    void Copy(const List<T>& c)
    {
      p.Size = 0;

      int i;
      for(i=0; i<c.p.Size; i++)
	{
	  this->Add(c.p.List[i]);
	}
    }

  protected:

    virtual bool IsValid(const T& w) const = 0;

  private:

    struct Pars
    {
      int Size;
      int BufferSize;
      T *List;
    }
    p;

  private:

    List(const List &);           // Not implemented
    void operator=(const List&);  // Not implemented.
  };
};

#endif  // __NG_LIST_H__
