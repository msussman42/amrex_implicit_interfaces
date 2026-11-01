#include <iostream>
struct Expression {
   float value;
   virtual void evaluate() = 0;
   virtual void derive(float seed) = 0;
};
struct Variable: public Expression {
   float partial;
   Variable(float value) {
      this->value = value;
      partial = 0.0f;
   }
   void evaluate() {}
   void derive(float seed) {
      partial += seed;
   }
};
struct Plus: public Expression {
   Expression &a, &b;
   Plus(Expression &a, Expression &b): a(a), b(b) {}
   void evaluate() {
      a.evaluate();
      b.evaluate();
      value = a.value + b.value;
   }
   void derive(float seed) {
      a.derive(seed);
      b.derive(seed);
   }
};
struct Multiply: public Expression {
   Expression &a, &b;
   Multiply(Expression &a, Expression &b): a(a), b(b) {}
   void evaluate() {
      a.evaluate();
      b.evaluate();
      value = a.value * b.value;
   }
   void derive(float seed) {
      a.derive(b.value * seed);
      b.derive(a.value * seed);
   }
};
int main () {
   // Example: Finding the partials of z = x * (x + y) + y * y at (x, y) = (2, 3)
   Variable x(2), y(3);
   Plus p1(x, y); Multiply m1(x, p1); Multiply m2(y, y); Plus z(m1, m2);
   z.evaluate();
   std::cout << "z = " << z.value << std::endl;
   // Output: z = 19
   z.derive(1);
   std::cout << "∂z/∂x = " << x.partial << ", "
             << "∂z/∂y = " << y.partial << std::endl;
   // Output: ∂z/∂x = 7, ∂z/∂y = 8
   return 0;
}
