namespace Algorhi{
  class Func{
    public:
      Func(int);
      Func(int a ,int b);
      Func(const Func& other);
      virtual ~Func();
      void get();

    private:
      int swap;
      int m,n;
  };
}
