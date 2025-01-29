clear all

% Parameters
numSamples = 1000;  % Total samples
numInputs = 2;      % Number of inputs
numOutputs = 2;     % Number of outputs
sequenceLength = 10; % Length of each sequence

% Generate sample time-series data
t = linspace(0, 10, numSamples)'; % Time vector
inputs = [sin(2*pi*0.1*t), cos(2*pi*0.1*t)]; % Two sine-wave inputs
outputs = [0.5 * inputs(:, 1) + 0.3 * inputs(:, 2) + 0.1 * randn(numSamples, 1), ...
           0.2 * inputs(:, 1).^2 - 0.1 * inputs(:, 2) + 0.05 * randn(numSamples, 1)];

% Normalize inputs and outputs for neural network training
inputs = normalize(inputs);
outputs = normalize(outputs);

% Reshape data into sequences
numSequences = floor(numSamples / sequenceLength);

% Reshape data to [sequenceLength, numSequences, numInputs] for inputs
inputs = reshape(inputs(1:numSequences*sequenceLength, :), sequenceLength, numSequences, numInputs);

% Reshape data to [sequenceLength, numSequences, numOutputs] for outputs
outputs = reshape(outputs(1:numSequences*sequenceLength, :), sequenceLength, numSequences, numOutputs);

% Convert data to cell arrays
trainInputs = squeeze(mat2cell(inputs, sequenceLength, ones(1, numSequences), numInputs))';
trainOutputs = squeeze(mat2cell(outputs, sequenceLength, ones(1, numSequences), numOutputs))';

% Verify the dimensions
disp(size(trainInputs)); % Should print [100, 1]
disp(size(trainInputs{1})); % Should print [10, 2]

% Split data into training (70%), validation (15%), and testing (15%)
numTrain = round(0.7 * numSequences);
numVal = round(0.15 * numSequences);

trainInputsFinal = trainInputs(1:numTrain);
trainOutputsFinal = trainOutputs(1:numTrain);

valInputs = trainInputs(numTrain+1:numTrain+numVal);
valOutputs = trainOutputs(numTrain+1:numTrain+numVal);

testInputs = trainInputs(numTrain+numVal+1:end);
testOutputs = trainOutputs(numTrain+numVal+1:end);

% Define the RNN layers
layers = [
    sequenceInputLayer(numInputs)                  % Input layer for sequences
    lstmLayer(50, 'OutputMode', 'sequence')        % LSTM layer with 50 hidden units
    fullyConnectedLayer(numOutputs)                % Fully connected layer for outputs
    regressionLayer];                              % Regression layer for system identification

% Training options
options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 16, ...
    'ValidationData', {valInputs, valOutputs}, ...
    'ValidationFrequency', 10, ...
    'InitialLearnRate', 0.001, ...
    'Plots', 'training-progress', ...
    'Verbose', true);

% Train the network
net = trainNetwork(trainInputsFinal, trainOutputsFinal, layers, options);

% Predict outputs for the test set
predictedOutputs = predict(net, testInputs);

% Convert cell arrays to matrices for comparison
actualOutputs = cell2mat(testOutputs);
predictedOutputsMat = cell2mat(predictedOutputs);

% Compute RMSE
rmse = sqrt(mean((actualOutputs - predictedOutputsMat).^2, 'all'));
fprintf('Test RMSE: %.4f\n', rmse);

% Plot comparison of actual and predicted outputs
figure;
subplot(2, 1, 1);
plot(actualOutputs(:, 1), 'b'); hold on;
plot(predictedOutputsMat(:, 1), 'r--');
title('Output 1: Actual vs Predicted');
legend('Actual', 'Predicted');

subplot(2, 1, 2);
plot(actualOutputs(:, 2), 'b'); hold on;
plot(predictedOutputsMat(:, 2), 'r--');
title('Output 2: Actual vs Predicted');
legend('Actual', 'Predicted');

