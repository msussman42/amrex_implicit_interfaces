# Code written by Chris Caballero
import torch
import math
import time
from torch.utils.tensorboard import SummaryWriter

class Net:
    # We will be using the Mean Squared Error (MSE) loss function
    loss_function = torch.nn.MSELoss(reduction='sum')
    # These are our inputs and outputs
    # y : expected_outputs
    # x : input_set
    input_set = torch.linspace(-math.pi, math.pi, 2000)
    outputs_expected = torch.sin(input_set)
    # Constructor
    def __init__(self):
        # y is a linear function of (x, x^2, x^3)
        # We prepare the tensor (x, x^2, x^3) : {x^p : p = 1, 2, 3}  
        p = torch.tensor([0,1,2,3])
        self.inputs = self.input_set.unsqueeze(-1).pow(p)

        # nn.Sequential : applies sub-modules in sequence
        # nn.Linear     : computes output using a linear function
        # nn.Linear     : holds internal tensors for its weight and bias
        # nn.Flatten    : flattens the output of the linear layer to a 1D tensor
        self.model = torch.nn.Sequential(
            torch.nn.Linear(4, 2),
            torch.nn.Linear(2, 1),
            torch.nn.Flatten(0, 1),
        )
        self.evaluate_model()
        

        # This is our predicted outputs
        self.outputs_predicted = self.model(self.inputs)

        # Determines the number of weights for the model and 
        # Copy the weights to keep track of them and perform operations
        self.num_weights = self.get_num_weights()
        for param in self.model.parameters():
            print(param)


        self.weights = []
        self.export_weights(self.weights)

        zero = []
        self.export_zeros(zero)
        self.zero_grad = self.evaluate_model(zero)[1]

    # Accessor Functions:
    # Counts the number of weights in the model
    def get_num_weights(self):
        num_weights = 0
        for param in self.model.parameters():
            num_weights += 1
        return num_weights
    
    # Copies the weights from the model into a given list
    def export_weights(self, weights):
        weights.clear()
        for param in self.model.parameters():
            weights.append(torch.flatten(param.clone()))
    
    # Copies weights from a given list into the model
    def import_weights(self, weights):
        i = 0
        with torch.no_grad():
            for param in self.model.parameters():
                param[0] = weights[i].clone()
                i += 1
    
    # Copies the gradients from the model into a given list
    def export_gradients(self, gradients):
        gradients.clear()
        for param in self.model.parameters():
            gradients.append(torch.flatten(param.grad.clone()))
    
    # Copies gradients from a given list into the model
    def import_gradients(self, gradients):
        i = 0
        with torch.no_grad():
            for param in self.model.parameters():
                param.grad[0] = gradients[i].clone()
                i += 1

    # First sets each element in the weights to 0
    # Then copies weights from the model into a given list
    def export_zeros(self, weights):
        weights.clear()
        for param in self.model.parameters():
            weights.append(0*torch.flatten(param.clone()))

    # Mutator Functions:
    # Updates the weights of the object to a new set of weights
    def update_weights(self, new_weights):
        self.weights.clear()
        i = 0
        for i in range(self.num_weights):
            self.weights.append(torch.flatten(new_weights[i].clone()))
        self.outputs_predicted = self.model(self.inputs)
    
    # Evaluates the model with a given set of weights
    # Defaults to current weights in the model
    def evaluate_model(self, weights = []):
        # evaluating the model with new weights, 
        # changes its internal weights and gradients,
        # so we store the originals to reset the model
        if weights:
            model_weights = []
            model_gradients = []
            self.export_weights(model_weights)
            self.export_gradients(model_gradients)
            self.import_weights(weights)

        # prediction : inputs evaluated using the current model state
        prediction = self.model(self.inputs)

        # loss struct contains the loss function evaluation and gradients
        loss_struct = self.loss_function(prediction, self.outputs_expected)
        loss_value = loss_struct.item()

        # backwards step computes the gradients of the model's weights
        self.model.zero_grad()
        loss_struct.backward()

        # store gradients for output
        loss_gradient = []
        self.export_gradients(loss_gradient)
        
        # if new weights were used to evaluate the model,
        # we reset the weights and gradients to their original state
        if weights:
            self.import_weights(model_weights)
            self.import_gradients(model_gradients)
        
        return loss_value, loss_gradient
    
    # Returns the loss using the predicted weights in the network
    def get_loss(self):
        return self.loss_function(self.outputs_predicted, self.outputs_expected)

    # Optimizes the weights of the model using:
    # Nonlinear Conjugate Gradient with Fletcher-Reeves
    def train(self, writer = None):
        print('Training using Nonlinear Conjugate Gradient:')
        # uses fletcher-reeves method to compute beta
        def compute_beta(delta, remainder):
            delta_old = delta
            delta = 0
            for i in range(self.num_weights):
                delta += torch.dot(remainder[i], remainder[i])
            return delta / delta_old, delta

        # loss : loss function value
        # remainder : remainder term (starts as gradient)
        # direction : direction of descent
        loss, remainder = self.evaluate_model()
        direction = []

        # delta is used for fletcher-reeves beta computation
        # delta = dot(remainder, remainder)
        delta = 0
        for i in range(self.num_weights):
            delta += torch.dot(remainder[i], remainder[i])

        for i in range(self.num_weights):
            # set remainder to -gradient
            remainder[i] = torch.neg(remainder[i])
            # set direction to -gradient (remainder)
            direction.append(remainder[i].clone())
        
        refresh_tracker = 0
        refresh_threshold = 10
        # tolerance = 1e-2 * len(self.input_set)
        tolerance = 9
        epoch = 0
        epochs = 5
        for epoch in range(epochs):
            print(epoch)
            # take a step in the rigth direction
            new_weights = self.line_search(remainder, direction)
            self.update_weights(new_weights)
            self.import_weights(self.weights)

            loss, remainder = self.evaluate_model()

            if loss <= tolerance:
                break

            for i in range(self.num_weights):
                remainder[i] = torch.neg(remainder[i])  
            # computes our beta and stores the updated delta
            beta_fr, delta = compute_beta(delta, remainder)

            # update direction for orthogonal search directions
            # descent monitors if our direction is a descent direction
            descent = 0
            for i in range(self.num_weights):
                direction[i] = remainder[i].add(beta_fr * direction[i])
                descent += torch.dot(remainder[i], direction[i])
 
            # reset direction = -gradient if:
            # threshold number of iterations completed, or
            # descent <= 0 (not going in descent direction)
            refresh_tracker += 1
            if refresh_tracker == refresh_threshold or descent <= 0:
                for i in range(self.num_weights):
                    direction[i] = remainder[i].clone()
                    refresh_tracker = 0
            epoch += 1
        #  print('\tTotal training iterations: {0}'.format(epoch))
        self.outputs_predicted = self.model(self.inputs)

    # Optimizes the weights using gradient descent
    def train_GD(self):
        print('Training using Gradient Descent:')
        epoch = 0
        epochs  = 5000
        learning_rate = 1e-6
        # tolerance = 1e-2 * len(self.input_set)
        tolerance = 19
        for epoch in range(epochs):
            loss_value = self.evaluate_model()[0]
            if loss_value < tolerance:
                break
            if epoch % 100 == 0:
                print(epoch, loss_value)

            with torch.no_grad():
                for param in self.model.parameters():
                    param -= learning_rate * param.grad
        print('\tTotal training iterations: {0}'.format(epoch))
        self.export_weights(self.weights)
        self.outputs_predicted = self.model(self.inputs)

    # Utility Functions:
    # Computes optimized step-size (alpha) in a given direction
    def line_search(self, remainder, direction):
        # offsets weights in a given direction with step-size alpha
        # offset weights : w + alpha * d
        def offset_weights(weights, alpha):
            for i in range(self.num_weights):
                weights[i] += alpha * direction[i].clone()
        
        # wolfe-condition (1) : f(offset weights) <= f(weights) + c1 * alpha * dot(f'(weights), direction)
        # wolfe-condition (2) : abs(dot(f(offset weights), direction)) <= -c2 * dot(f'(weights), direction)
        # out(boolean) : returns true if both conditions are met, false otherwise
        def wolfe_conditions(candidate, alpha, c1 = 1e-4, c2 = 1e-2):
            # evaluates the model using the offset weights
            # stores the respective loss value and gradients 
            # f(offset weights), f'(offset weights)
            loss_candidate, loss_candidate_grad = self.evaluate_model(candidate)
            # evaluates the model using its current weights
            # stores the respective loss value and gradients 
            # f(weights), f'(weights)
            loss, loss_grad = self.evaluate_model()

            # left1 : f(offset weights)
            # right1 : f(weights)
            left1 = loss_candidate
            right1 = loss

            # print(left1, right1)
            # print(loss_candidate_grad, loss_grad)

            for i in range(self.num_weights):
                right1 += c1 * alpha * torch.dot(loss_grad[i], direction[i]).item()

            return left1 <= right1

        # Better Alpha?
            # Try to use 'optimal' alpha initially 
            # I set gamma = dot(remainder, remainder)/dot(direction, mul(A, direction))
            # For certain systems, this converges in n steps, very fast
            # It seems inefficient to have to run the model for this (if the system does not converge in n steps)
        direction_grad = self.evaluate_model(direction)[1]
        A_direction = []
        for i in range(self.num_weights):
            A_direction.append(torch.add(direction_grad[i], self.zero_grad[i]))

        num = 0
        denom = 0
        for i in range(self.num_weights):
            num += torch.dot(remainder[i], remainder[i]).item()
            denom += torch.dot(direction[i], A_direction[i]).item()
        gamma = num / denom
        alpha = gamma ** 0
        #     # small starting step size (maybe better in other cases?)
        # gamma = 5e-1
        # alpha = gamma**0

        # candidate list : stores the offset weights
        # we start the search using the current weights 
        candidate = []
        self.export_weights(candidate)
        conditions_satisfied = wolfe_conditions(candidate, alpha)
        # offset_weights(candidate, alpha)
        # if gamma >= 1:
        #     gamma = 5e-1
        
        p = 1
        max_search_iterations = 20
        while not conditions_satisfied and p < max_search_iterations:
            alpha = gamma ** p
            self.export_weights(candidate)
            offset_weights(candidate, alpha)
            conditions_satisfied = wolfe_conditions(candidate, alpha)
            print(candidate)
            p += 1
        return candidate
            
    # Prints the weights in the model
    def show_model_weights(self):
        print('The weights in the model:')
        for param in self.model.parameters():
            print('\t', param[0])
    
    # Prints the weights in the class
    def show_object_weights(self):
        print('The weights in the object:')
        for i in range(self.num_weights):
            print('\t', self.weights[i])
    
    # Prints the gradient of the weights in the model
    def show_model_gradients(self):
        print('The gradients in the model:')
        for param in self.model.parameters():
            print('\t', param.grad)

# Prints the weights given as a parameter
def show_weights(weights):
    print('The weights are:')
    for i in range(len(weights)):
        print('\t', weights[i])

# Prints the gradients given as a parameter
def show_gradients(gradients):
    print('The gradients are:')
    for i in range(len(gradients)):
        print('\t', gradients[i])


def main():   
    # net = Net()
    # print('\tInitial Loss Value: {0}'.format(net.get_loss()))
    # net.train()
    # print('\tFinal Loss Value: {0}'.format(net.get_loss()))

    comp = Net()
    print('\tInitial Loss Value: {0}'.format(comp.get_loss()))
    comp.train_GD()
    print('\tFinal Loss Value: {0}'.format(comp.get_loss()))


if __name__ == "__main__":
    main()

# SANITY CHECKS:
    # IMPORT EXPORT SANITY CHECK:
        ### WEIGHTS
            # # import model weights to w
            # net.import_weights(w)
            # # show 
            # net.show_object_weights()
            # net.show_model_weights()

            # w[0] -= 2*w[0]

            # net.show_object_weights()
            # net.show_model_weights()
            # # export w to model
            # net.export_weights(w)
            # net.show_model_weights()

            # net.update_weights(w)
            # net.show_object_weights()

        ### GRADIENTS
            # l, g = net.evaluate_model()
            # net.show_model_gradients()
            # show_gradients(g)

            # l2, g2 = net.evaluate_model(w)
            # net.show_model_gradients()
            # show_gradients(g2)

            # net.import_gradients(g2)
            # net.show_model_gradients() 
    
    # COMPARISONS
        # comp = Net()
        # start_GD = time.perf_counter()
        # comp.train_GD()
        # end_GD = time.perf_counter()
        # comp.show_object_weights()

        # print('Time to train using NCG: {0} seconds'.format(end - start))
        # print('Time to train using GD: {0} seconds'.format(end_GD - start_GD))

    # MODEL EVALUATION SANITY CHECK:
        # show_weights(w)
        # w[0] -= 2*w[0]

        # loss, grad = net.evaluate_model(w)
        # print(loss, grad)
        
        # show_weights(w)
        # net.show_model_weights() 
        # net.show_model_gradients()   
