import numpy as np
import matplotlib.pyplot as plt
from sklearn import tree


class Waveform(object):
    def __init__(self,windspeed,waveform):
        self.windspeed = windspeed
        self.waveform = waveform
        self.rms = np.sqrt(np.mean(np.square(waveform)))




class WaveformAnalyzer(object):
    """Provides a framework of tools for analyzing the Willoughby waveforms

    Attributes
     data: waveforms in np format
    
    """
    def __init__(self):
        """load the willoughby file"""

        print "Loading massive Willoughby dataset"

        #load the data,
        #self.data = np.loadtxt('will56530.txt', usecols=range(9,608))
        #self.data = np.loadtxt('will56530.txt')
        
        #for w in range(len(self.data)-1):
            



        #remove the first 9 elements of every waveform array 
        #to get rid of unnessary info
        #for w in range(len(self.data)-1):
         #   self.waveforms[w][:] = self.data[w][9:len(self.data[w])]
            
            #self.waveforms[w] = wave
            #print self.waveforms[0]
            
            

    def plot_waveform(self,waveform_ids):
        """Take a specified waveform and plot it """

        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])
        
        #for x in range(len(waveform_ids)-1):
        #plt.plot(self.data[waveform_ids[x]][6:])

        plt.ylim(0,0.2)
        plt.plot(self.data[0][6:])
        #plt.plot(self.data[29][6:])
        #plt.plot(self.data[500][6:])
        plt.plot(self.data[56660][6:])
            
        plt.ylabel('Waveform intensity')
        plt.xlabel('Index')        
        plt.show()
        
if __name__ == '__main__':


    print "Loading data"
    data = np.loadtxt('will56530.txt')

    waveform_database = []
    
    print "Formatting data"

    for d in data:
        w = Waveform(d[3], np.array(d[8:608]))
        waveform_database.append(w)


    waveform_database = np.array(waveform_database)


    #Take half of the Willoughby storms and 
    #use them as training data
    np.random.shuffle(waveform_database)

    half_index = int((len(waveform_database)/2))
    train = waveform_database[1:half_index]
    test = waveform_database[half_index+2:]

    X = [] #training samples
    y = [] #training labels

    #Get the waveforms into a training sample array
    for t in train:
        X.append(t.waveform)
        
    for u in test:
        y.append(u.windspeed)
    
    X = np.array(X)
    y = np.array(y)

    print X.shape
    print X[0][0]

    clf = tree.DecisionTreeRegressor()     

    print "Training model"
    clf = clf.fit(X,y)


    final_predictions = []
    perror = []
    

    #Run some tests
    for x in test:
        prediction = clf.predict(x.waveform)[0]
        final_predictions.append(prediction)
        actual = x.windspeed

        pe = abs(prediction-actual)/actual#percentage error
        
        perror.append(pe)

        print prediction, actual, pe

        

     
    final_predictions = np.array(final_predictions)
    final_predictions = np.reshape(final_predictions, -1)
    print final_predictions[1:10]
    print final_predictions.shape

    perror = np.array(perror)
    mean_error = np.mean(perror)

    print mean_error
    

    
                 
             

        

