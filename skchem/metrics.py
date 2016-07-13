#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

def bedroc_score(y_true, y_pred, decreasing=True, alpha=20.0):

    """BEDROC metric implemented according to Truchon and Bayley.

    The Boltzmann Enhanced Descrimination of the Receiver Operator
    Characteristic (BEDROC) score is a modification of the Receiver Operator
    Characteristic (ROC) score that allows for a factor of *early recognition*.

    References:
        The original paper by Truchon et al. is located at `10.1021/ci600426e
        <http://dx.doi.org/10.1021/ci600426e>`_.

    Args:
        y_true (array_like):
            Binary class labels. 1 for positive class, 0 otherwise.
        y_pred (array_like):
            Prediction values.
        decreasing (bool):
            True if high values of ``y_pred`` correlates to positive class.
        alpha (float):
            Early recognition parameter.

    Returns:
        float:
            Value in interval [0, 1] indicating degree to which the predictive
            technique employed detects (early) the positive class.
     """

    assert len(y_true) == len(y_pred), \
     'The number of scores must be equal to the number of labels'

    N = len(y_true)
    n = sum(y_true == 1)

    if decreasing:
        order = np.argsort(-y_pred)
    else:
        order = np.argsort(y_pred)

    m_rank = (y_true[order] == 1).nonzero()[0]

    s = np.sum(np.exp(-alpha * m_rank / N))

    r_a = n / N

    rand_sum = r_a * (1 - np.exp(-alpha))/(np.exp(alpha/N) - 1)

    fac = r_a * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha/2 - alpha * r_a))

    cte = 1 / (1 - np.exp(alpha * (1 - r_a)))

    return s * fac / rand_sum + cte
