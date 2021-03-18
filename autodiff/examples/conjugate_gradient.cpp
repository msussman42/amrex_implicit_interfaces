#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


//Given Functions
	//J(x)
	//grad_J(x)
	//transpose(V)
	//max()
	//dotproduct(V,V)
	//mult(V,V)
	//Show(V)
void line_search(vector<double> x, vector<double> p, double &alpha);

int main() {

// Create Tensors to hold input and outputs.
x = torch.linspace(-math.pi, math.pi, 2000)
y_supervised = torch.sin(x)  # supervised learning
p = torch.tensor([1, 2, 3])
xx = x.unsqueeze(-1).pow(p)

model = torch.nn.Sequential(
    torch.nn.Linear(3, 1),
    torch.nn.Flatten(0, 1)
)

loss_fn = torch.nn.MSELoss(reduction='sum')

ncomp_weights=0
get_number_of_weights(ncomp_weights)
allocate_array(xCG,ncomp_weights)
copy_weights_from_model(xCG)
allocate_array(xCG_output,ncomp_weights)
call optimize_weights_using_FletcherReeves(xCG,ncomp_weights,
		y_supervised,xx,xCG_output);

}


optimize_weights_using_FletcherReeves(
		input: xCG,input: ncomp_weights,
		input: y_supervised,
		input: xx,
		output: xCG_output) {

 int n = ncomp_weights;
 //NON-LINEAR CG
 //J(x,true)
 // double epsilon = pow(10, -5);

 //alpha (what does alpha do)
 double a;
 //direction vector
 vector<double> d;
 //current x_i
 vector<double> x(0,n);
 for (int i=0;i<n;i++) {
	 x[i]=xCG[i];
 }
 //Fletcher-Reeves Beta for determining conjugate directions
 vector<double> B_fr;
 //current and previous remainder terms
 vector<double> r;
 vector<double> r_1;
 vector<double> gradJ; // n components

 model_eval_and_grad(y_supervised,x,xx,loss_function_struct,
   gradJ);

 r = -gradJ;
 d = r;
 b = r

//get loss and gradient
//on pytorch

 int j = 0;
 while(norm_2(gradJ) > tolerance) {

  //find a[i] to minimize f(x[i] + a[i]d[i])
  //not sure (should do line search???)
  //a = -dotproduct(transpose(r), r)/dotproduct(transpose(d), (b-grad_J(d)));
  line_search(y_supervised,x,xx,d,a);

  //find next x_i and next remainder term w/ new x
  x = x + mult(a, d);
  model_eval_and_grad(y_supervised,x,xx,loss_function_struct,
   gradJ);
  r_1 = -gradJ;

  //fletcher-reeves
  B_fr = dotproduct(transpose(r_1), r_1)/dotproduct(transpose(r), r);
  d = r_1 + mult(B_fr, d);

  r = r_1, j++;
 }
 xCG_output=x;

Show(r);
//Stochastic CG
//J(x,false)
int T = 30;
int m = 10;
//given w0;
vector<double> d0;
vector<double> g;
vector<double> g1;
vector<double> h;
vector<double> w;
vector<double> u;
vector<double> x0;
vector<double> x1;
vector<double> S;



h = grad_J(w);
for(int k = 0; k < T-1; k++) {
	u = grad_J(w);
	x1 = x0 = w;
	g = h;
	d0 = -g;

	for(int t = 0; t < m-1; t++) {
		double a;
		//Sample a mini-batch
		//S = new vector<double>(Random() % n, 0);
		
		//figure out how to do this^^^
			// for(int i = 0; i < S.size(); i++) {
			// 	//how do i select this properly??
			// 	S[i] = Random() % n;
			// }
		line_search(x,d,a);
		x1 += a*d;
		//not sure how else to get gradient 
		auto g0 = g1;
		g1 = grad_J(x1) - grad_J(x0) + u;
		//Fletcher-Reeves
		B_fr = dotproduct(transpose(g1), g1)/dotproduct(transpose(g0), g0);

		d = -g1 + B_fr*d;
	}
	h = g1;
	w = x1;
}

	return 0;
}

void line_search(y_supervised,vector<double> x,
	       xx,vector<double> d, double &alpha) {

	alpha = (double) rand() % MAX_RAND;


	double step_size;
	double c = (double) rand() % MAX_RAND;
	vector<double> x_candidate=x+pow(alpha,i)*d;
        call model_eval_and_grad(y_supervised,x_candidate,
			xx,loss_function_struct,
                        gradJ);
        call model_eval_and_grad(y_supervised,x,
			xx,loss_function_struct_base,
                        gradJ_base);
	double left = loss_function_struct.item();
	double right = loss_function_struct_base.item() + 
		c*pow(alpha,i)*dotproduct(transpose(d), gradJ_base);


	int i = 0, max_itr = 1000;
	while (i < max_itr && left > right) {
		i++;


 	 vector<double> x_candidate=x+pow(alpha,i)*d;
         call model_eval_and_grad(y_supervised,x_candidate,
			xx,loss_function_struct,
                        gradJ);
         call model_eval_and_grad(y_supervised,x,
			xx,loss_function_struct_base,
                        gradJ_base);
	 double left = loss_function_struct.item();
	 double right = loss_function_struct_base.item() + 
		c*pow(alpha,i)*dotproduct(transpose(d), gradJ_base);

	}
	if(left <= right) {
		alpha = pow(alpha,i);
	}
}


/*
	//d0 = -gradient0
	//ak = -gradientk*dk/(transpose(dk)*Q*dk)
	//xk+1 = xk + ak*dk
	//gk = Q*xk - b
	//dk+1 = -gradientk+1 + Bk*dk
	//Bk = trans(gradientk+1)*Q*dk/(transpose(dk)*Q*dk)
	// int i = 0;
	// //LINEAR CG
	// for(i; i < n; i++) {
	// 	//not sure about this (J for A????)
	// 	a[i] = dotproduct(transpose(r[i]),r[i])/dotproduct(transpose(d[i]),J(d[i])+b);
	// 	x[i+1] = x[i] + a[i]*d[i];
	// 	r[i+1] = r[i] - a[i]*J(d[i]);
	// 	if(r[i+1] < epsilon) {
	// 		break;
	// 	}
	// 	//fletcher reeves
	// 	B[i+1] = dotproduct(transpose(r[i+1]),r[i+1])/dotproduct(transpose(r[i]),r[i]);

	// 	d[i+1] = r[i+1] + B[i+1]*d[i];
	// }

		//polak-ribiere
		//B[i+1] = max(transpose(r[i+1])*(r[i+1]-r[i])/(transpose(r[i])*r[i]), 0);

		// d[i+1] = r[i+1] + B[i+1]*d[i];

		// if(++k == m || dotproduct(transpose(r[i]),d[i]) <= 0) {
		// 	d[i+1] = r[i+1];
		// 	k = 0;
		// }
