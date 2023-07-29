#!/Softwares/anaconda3/envs/md/bin/python


                                                    #**********************************************************#
                                                    #                                                          #
                                                    #                 OWNER : SUBINOY ADHIKARI                 #
                                                    #                                                          #
                                                    #**********************************************************#

import sys
import os
import pickle
import tensorflow as tf
from tensorflow.keras import Model
from tensorflow.keras.layers import Input, Dense, BatchNormalization, Activation
from tensorflow.keras.losses import MeanSquaredError
from tensorflow.keras.optimizers import Adam
import numpy as np

def gpu_id(gpu):
	os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
	os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu) # Model will be trained on the specified GPU #




class dense_autoencoder:
    """
    This is a fully connected multi layer perceptron.
    """

    def __init__(self,
                 input_shape,
                 encoder_neurons,
                 encoder_activation,
                 encoder_wt_initialization,
                 latent_neurons,
                 latent_activation,
                 latent_wt_initialization,
                 decoder_neurons,
                 decoder_activation,
                 decoder_wt_initialization):

        self.input_shape=input_shape
        self.encoder_neurons=encoder_neurons
        self.encoder_activation=encoder_activation
        self.encoder_wt_initialization=encoder_wt_initialization
        self.latent_neurons=latent_neurons
        self.latent_activation=latent_activation
        self.latent_wt_initialization=latent_wt_initialization
        self.decoder_neurons=decoder_neurons
        self.decoder_activation=decoder_activation
        self.decoder_wt_initialization=decoder_wt_initialization

        # Number of layers in the encoder and decoder
        self._num_encoder_layers=len(encoder_neurons)
        self._num_decoder_layers=len(decoder_neurons)


        self.encoder=None
        self.decoder=None
        self.autoencoder=None
        self.encoder_input=None

        # Call the dense autoencoder model
        self._dense_model()

    # The dense autoencoder model
    def _dense_model(self):
        self._encoder_model()
        self._decoder_model()
        self._autoencoder_model()

    # --------------------------------ENCODER--------------------------------#

    def _encoder_model(self):
        enc_inp=self._encoder_input() # Input layer of the encoder, with neurons equal to the number of input features #
        enc_lyr=self._encoder_layers(enc_inp) # Stacking of the hidden layers of the encoder, excluding the latent layer #
        lec=self._latent_layer(enc_lyr) # Latent layer #
        self.encoder=Model(inputs=enc_inp,
                           outputs=lec,
                           name="ENCODER") # Encoder #
        self.encoder_input=enc_inp

    # Encoder input layer #
    def _encoder_input(self):
        enc_inp=Input(shape=self.input_shape,
                      name="ENCODER_INPUT")
        return enc_inp

    # Stacking the hidden layers of the encoder #
    def _encoder_layers(self, enc_inp):
        enc = enc_inp
        enc=BatchNormalization(name="BATCH_NORMALIZATION_FOR_THE_INPUT_LAYER_OF_THE_ENCODER")(enc)
        for index in range(self._num_encoder_layers):
            enc=Dense(units=self.encoder_neurons[index],
                      kernel_initializer=self.encoder_wt_initialization[index],
                      use_bias=False,
                      name=f"ENCODER_LAYER_{index+1}")(enc)

            enc=BatchNormalization(name=f"BATCH_NORMALIZATION_FOR_ENCODER_LAYER_{index+1}")(enc)
            enc=Activation(self.encoder_activation[index])(enc)
        return enc

    # Latent layer #
    def _latent_layer(self, enc_lyr):
        lec=Dense(units=self.latent_neurons,
                  kernel_initializer=self.latent_wt_initialization,
                  use_bias=False,
                  name="LATENT_LAYER")(enc_lyr)

        lec=BatchNormalization(name=f"BATCH_NORMALIZATION_FOR_LATENT_LAYER")(lec)
        lec=Activation(self.latent_activation)(lec)
        return lec

    # --------------------------------DECODER--------------------------------#

    def _decoder_model(self):
        dec_inp=self._decoder_input() # Input layer of the decoder, with neurons equal to the number of latent neurons #
        dec_lyr=self._decoder_layers(dec_inp) # Stacking of the hidden layers of the decoder, excluding the output layer #
        dec_out=self._output_layer(dec_lyr) # Output layer #
        self.decoder=Model(inputs=dec_inp,
                           outputs=dec_out,
                           name="DECODER") # Decoder #

    # Decoder input layer #
    def _decoder_input(self):
        dec_inp=Input(shape=self.latent_neurons,
                      name="DECODER_INPUT")
        return dec_inp

    # Stacking the hidden layers of the decoder #
    def _decoder_layers(self, dec_inp):
        dec=dec_inp
        dec=BatchNormalization(name="BATCH_NORMALIZATION_FOR_THE_INPUT_LAYER_OF_THE_DECODER")(dec)
        for index in range(self._num_decoder_layers):
            dec=Dense(units=self.decoder_neurons[index],
                      kernel_initializer=self.decoder_wt_initialization[index],
                      use_bias=False,
                      name=f"DECODER_LAYER_{index+1}")(dec)

            dec=BatchNormalization(name=f"BATCH_NORMALIZATION_FOR_DECODER_LAYER_{index+1}")(dec)
            dec=Activation(self.decoder_activation[index])(dec)
        return dec

    # Output layer #
    def _output_layer(self, dec_lyr):
        dec_out=Dense(units=self.input_shape,
                      kernel_initializer='glorot_uniform',
                      use_bias=True,
                      name="DECODER_OUTPUT")(dec_lyr)
        dec_out=Activation('sigmoid')(dec_out)
        return dec_out

    # --------------------------------AUTOENCODER--------------------------------#
    def _autoencoder_model(self):
        input_model=self.encoder_input
        output_model=self.decoder(self.encoder(input_model))
        self.autoencoder=Model(inputs = input_model,
                               outputs = output_model,
                               name = "AUTOENCODER")

    # --------------------------------SUMMARY OF THE AUTOENCODER--------------------------------#
    def summary(self):
        self.encoder.summary()
        self.decoder.summary()
        self.autoencoder.summary()
    # --------------------------------COMPILING THE AUTOENCODER--------------------------------#
    def compile(self,
                optimizer=Adam,
                learning_rate=0.001,
                loss=MeanSquaredError()):
        self.autoencoder.compile(optimizer=optimizer(learning_rate),
                                 loss=loss)

# --------------------------------TRAINING THE AUTOENCODER--------------------------------#
    def train(self,
              x_train=None,
              x_test=None,
              y_train=None,
              y_test=None,
              batch_size=100,
              epochs=100,
              shuffle=True):
        self.autoencoder.fit(x=x_train,
                             y=x_train,
                             batch_size=batch_size,
                             epochs=epochs,
                             shuffle=shuffle,
                             validation_data=(x_test, x_test)) # x_train=training data, x_test=testing data, y_train=labels of the training data, y_test=labels of the testing data #

# --------------------------------SAVING THE AUTOENCODER--------------------------------#
    def save(self, save_dir_name="."):

        # Create the directory to save the autoencoder model #
        check_dir = os.path.isdir(save_dir_name) # Checks if the specified directory exists #
        if not check_dir: # If directory does not exist, then create it #
            os.makedirs(save_dir_name, exist_ok=False)
            print("Created directory : ", save_dir_name)
        else:
            os.makedirs(save_dir_name, exist_ok=True) # If directory exists, then overwrite it #
            print("Overwritten directory : ", save_dir_name)

       # Saving the parameters, weights and loss of the autoencoder model #
        parameters = [self.input_shape,
                      self.encoder_neurons,
                      self.encoder_activation,
                      self.encoder_wt_initialization,
                      self.latent_neurons,
                      self.latent_activation,
                      self.latent_wt_initialization,
                      self.decoder_neurons,
                      self.decoder_activation,
                      self.decoder_wt_initialization]
        parameters_filename=os.path.join(save_dir_name, "parameters.pkl") # Path and name of the file where the parameters would be saved #
        with open(parameters_filename, "wb") as file:
            pickle.dump(parameters, file)

        # Saving the weights of the autoencoder #
        weights_filename=os.path.join(save_dir_name, "weights.h5") # Path and name of the file where the weights would be saved #
        self.autoencoder.save_weights(weights_filename)

        # Saving the training loss and validation loss #
        training_loss_filename = os.path.join(save_dir_name, "training_loss.dat") # Path and name of the file where the training loss would be saved #
        np.savetxt(training_loss_filename, self.autoencoder.history.history['loss']) # Training loss #

        validation_loss_filename = os.path.join(save_dir_name, "validation_loss.dat") # Path and name of the file where the validation loss would be saved #
        np.savetxt(validation_loss_filename, self.autoencoder.history.history['val_loss']) # Validation loss #


    # --------------------------------LOADING THE AUTOENCODER--------------------------------#
    def _load_weights(self, weights_filename):
        self.autoencoder.load_weights(weights_filename)

    @classmethod
    def load(cls, save_dir_name="."):
        # Load the parameters #
        parameters_filename=os.path.join(save_dir_name, "parameters.pkl")
        with open(parameters_filename, "rb") as file:
            parameters=pickle.load(file)
        autoencoder=dense_autoencoder(*parameters)

        # Load the weights #
        weights_filename=os.path.join(save_dir_name, "weights.h5")
        autoencoder._load_weights(weights_filename)
        return autoencoder

    # --------------------------------RECONSTRUCTION--------------------------------#
    def reconstruct(self, data):
        latent_data=self.encoder.predict(data) # Get the latent data #
        reconstructed_data=self.decoder.predict(latent_data) # Get the reconstructed data #
        return latent_data, reconstructed_data
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # -----ENCODER-----#
    encoder_neurons = [512, 256, 128, 64, 32]
    encoder_activation = ['elu', 'elu', 'elu', 'elu', 'elu']
    encoder_wt_initialization = ['he_normal', 'he_normal', 'he_normal', 'he_normal', 'he_normal']

    # -----LATENT-----#
    latent_neurons = 2
    latent_activation = 'elu'
    latent_wt_initialization = 'he_normal'

    # -----DECODER-----#
    decoder_neurons = [32, 64, 128, 256, 512]
    decoder_activation = ['elu', 'elu', 'elu', 'elu', 'elu']
    decoder_wt_initialization = ['he_normal', 'he_normal', 'he_normal', 'he_normal', 'he_normal']

    # -----AUTOENCODER-----#
    ae = dense_autoencoder(input_shape=786,
                           encoder_neurons=encoder_neurons,
                           encoder_activation=encoder_activation,
                           encoder_wt_initialization=encoder_wt_initialization,
                           latent_neurons=latent_neurons,
                           latent_activation=latent_activation,
                           latent_wt_initialization=latent_wt_initialization,
                           decoder_neurons=decoder_neurons,
                           decoder_activation=decoder_activation,
                           decoder_wt_initialization=decoder_wt_initialization)
    ae.summary()
