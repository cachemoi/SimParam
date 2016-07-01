import numpy as np
import matplotlib.pyplot as plt

def MyLogRand (mu, sigma):

    sample = np.random.lognormal(mu, sigma, 1000)

    print(str(sample))

    count, bins, ignored = plt.hist(sample, 100, normed=True, align='mid')
    x = np.linspace(min(bins), max(bins), 10000)
    pdf = (np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2)) / (x * sigma * np.sqrt(2 * np.pi)))

    plt.plot(x, pdf, linewidth=2, color='r')
    plt.xlim(xmax=20)
    # plt.axis('tight')
    plt.show()

MyLogRand(1, 1)