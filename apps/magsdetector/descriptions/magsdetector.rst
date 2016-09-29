magsdetector detects events by finding similarities in defined master events
based on cross-correlation of waveforms or envelopes. An event (origin) is declared if
the waveform fit is above a configured threshold. Location and depth are taken
from the master events whereas time is taken from the time window currently
processed and magnitude is computed from the ratio of the master event PGV to
the continuous PGV.

.. math::
   :label: M

   M = M_{master} + \frac{1}{N} \sum_{i=1}^{N} \ln \frac{PGV_{continuous_i}}{PGV_{master_i}}

magsdetector can be ran in realtime or offline.

Workflow
========

magsdetector subscribes to all configured :confval:`channels`. On each received
record the common time window of all channels is computed and if sufficent data
are available they are processed.

.. figure:: media/magsdetector/stepping.png

   Processing of a data time window of contiuous channels. The step of each
   processing step is the time of one sample. See :ref:`mags-processing` for how
   each time window is processed.

Each processed time window is :ref:`matched <cross-correlation>` with the master
event. If the overall fit exceeds the :confval:`configured threshold <detector.threshold>` the detector
is started to search for the maximum fit within the
:confval:`configured timewindow <detector.window>`. The timestamp of the maximum
fit is taken as event time and an origin with a magnitude is sent to the
messaging.


.. _mags-processing:

Processing
==========

The data are processed according to the configuration. Both schemas are
illustrated below.

.. figure:: media/magsdetector/processing.png

   Processing schema in either frequency domain and time domain. **data** are
   the input data received from records and sensitivity corrected. **data** are
   also the outputs used by subsequent steps.


The processing is applied to the data of the master events as well as to the
continuous data. Processing in the frequency domain is usually more accurate and
does not introduce phase delays as with recursive filters. On the other hand
time domain processing is much faster.

Envelope
--------

The cross-correlation can be done on the seismograms itself or on
the envelopes of the seismograms. We compute the envelope using
the Hilbert-transform in frequency domain:

.. math::
   :label: Env

   \hat{y}_i^j = \sqrt{(y_i^j)^2 + (H\{y_i^j\})^2}

or a running root mean square in time domain:

.. math::
   :label: EnvRT

   \hat{y}_i^j = \sqrt{ \frac{2}{N} \sum_{k = i-N}^i (y_k^j)^2 }

Here y\ :sub:`i`\ :sup:`j` is the output signal of the single trace operation,
which in this case is the envelope of the signal. y\ :sub:`i` is the input signal
of the single trace operation, which in this case is the filtered seismogram.
i and j are counters for time and trace number, respectively. H(y)
is the Hilbert-transform of function y. N is the number of samples used for
the RMS which is calculated as:

.. math:: N = \frac{\text{sampling frequency}}{\text{envelope.hiFreq}}

Using the envelope of seismograms instead of the seismograms itself, the
cross-correlation becomes less sensitive to small changes in the source location
and in source mechanism between master signal and the earthquake to be detected.

Logarithm of trace
------------------

As an option we can compute the logarithm of a trace:

.. math::
   :label: Log

   y_i^j = sgn(y_i^j) \: ln|y_i^j|

Application of this single trace operation increases the importance
of small amplitudes in the cross-correlation. Therefore, it will gen-
erally amplify the noise level, which is an unwanted effect. On
the other hand using this option we can increase the importance
of small amplitude waves like coda waves in comparison to direct
P- and S-waves or we can increase the importance of recordings at
distant stations in comparison to close station (if we compute the
network cross-correlation).

Noise removal
-------------

If we use seismogram envelopes, the constant background noise
level causes a problem in the computation of the cross-correlation.
Even if we remove the mean value in the seismograms, a constant
offset appears in the envelopes. This offset value in the envelopes
corresponds to the standard deviation of noise in the original seismograms.
If we correlate a constant (noise) trace with the master
event, a large correlation coefficient results possibly causing wrong
detections. A simple solution to this problem is to remove the noise
level from the master event traces as well as from the current time
window traces. We estimate the noise level y\ :sup:`j` in a time window
j before the potential signal at envelope trace y\ :sub:`i`\ :sup:`j`:

.. math::
   :label: Noise

   \overline{y^j} = \frac{1}{N_n}\sqrt{\sum_{i=1}^{N_n} y_i^j}

Here i is the counter for time samples, N\ :sub:`n` is the number of time
samples in the time window to estimate the noise level, and j is a
counter for the traces (seismometers times components).
Then we remove the noise level, which is a constant offset in
the envelopes:

.. math::
   :label: NoiseRM

   \hat{y}_i^j = y_i^j - \overline{y^j}

Here, the time window, where we remove the offest (and where we
later compute the correlation) is generally different from the time
window, where we estimate the noise level.

.. _cross-correlation:

Cross-correlation
=================

Single trace
------------

In a first step we compute the single trace cross-correlation with
zero lag at channel j:

.. math::
   :label: Rj

   R^j = \frac{1}{\sqrt{(\sum_{i=1}^N x_i^j x_i^j) \: (\sum_{i=1}^N y_i^j y_i^j)}} \: \sum_{i=1}^N x_i^j y_i^j


Here i is the counter for time samples, N is the number of time samples in the
used correlation time window, x\ :sup:`j` is the filtered seismogram or envelope
of the master event after single trace operations described in :ref:`mags-processing`,
y\ :sup:`j` is the filtered seismogram or envelope of the current time window
after applying the same operations as to the master traces, and j is a counter
for the traces (seismometers times components).

This cross-correlation coefficient is computed for all traces j=1 to M, where
M/3 is the total number of seismometers of the master event, if we use
three-component sensors. The first condition for a detection is, that

.. math::
   :label: cond1

   R^j > \text{detector.channelThreshold, for all } j = 1 \text{ to } M_{min} \leq M

where R\ :sup:`j` is defined in equation :eq:`Rj`, R\ :sub:`1` is a configurable input
threshold to be chosen in the order of R\ :sub:`1` ≈ 0.6 - 0.8, and M\ :sub:`min`
is the minimum number of seismograms, where the earthquake is
detectable.

In other words: This condition requires that the current time window
and the master event shows similar waveforms (R\ :sup:`j` > R\ :sub:`1`) for at
least M\ :sub:`min` seismograms. We use the parameter M\ :sub:`min` to account for the
fact that some seismograms may show large local noise (R\ :sup:`j` < R\ :sub:`1`),
when a small earthquake occurs.

M\ :sub:`min` is also useful, when data transmission of some sensors is
interrupted. In that case the correlation coefficient is set to zero
R\ :sup:`j` = 0 < R\ :sub:`1` and the result is similar to a local noise disturbance.
A detection is only possible, if the number of noisy or otherwise
disturbed traces does not exceed M - M\ :sub:`min`.

M\ :sub:`min` is implicitely given through :confval:`detector.minimumChannelRatio`.

Network
-------

In a next step we compute the network cross-correlation with
zero lag, where two options, trace (A) and total (B), are available. What option
is used is configurable with :confval:`processing.normalization`.

.. note::

   The network cross-correlation is only computed if the ratio of time
   window traces and master event traces does not fall below :confval:`detector.minimumChannelRatio`
   and if the ratio of time window stations and master event stations does not fall
   below :confval:`detector.minimumStationRatio`. Otherwise the network cross-correlation
   is set to 0.

In option A we compute the average single trace cross-correlation coefficient
of the M\ :sub:`min` traces with the highest single trace cross-correlation
coefficient:

.. math::
   :label: RA

   R_{trace} &= \frac{1}{M_{min}} \sum_{j=1}^{M_{min}} R^j \\
       &= \frac{1}{M_{min}} \sum_{j=1}^{M_{min}} \frac{1}{\sqrt{(\sum_{i=1}^N x_i^j x_i^j) \: (\sum_{i=1}^N y_i^j y_i^j)}} \: \sum_{i=1}^N x_i^j y_i^j

In option B we compute the matrix cross-correlation with zero lag, where one
dimension of the matrix is time and the other dimension are the different traces:

.. math::
   :label: RB

   R_{total} = \frac{1}{\sqrt{(\sum_{j=1}^{M_{min}} \sum_{i=1}^N x_i^j x_i^j) (\sum_{j=1}^{M_{min}} \sum_{i=1}^N y_i^j y_i^j)}} \sum_{j=1}^{M_{min}} \sum_{i=1}^N x_i^j y_i^j

The major difference between R\ :sub:`trace` and R\ :sub:`total` is that in equation
:eq:`RA` the relative amplitudes between stations and components are neglected,
whereas they are taken into account in equation :eq:`RB`.

Then, the second and third condition for a detection are:

.. math::
   :label: cond2

   R_{trace} > \text{detector.threshold}

.. math::
   :label: cond3

   R_{total} > \text{detector.threshold}

Here :confval:`detector.threshold` is a configurable threshold, which judges the
network similarity and should be set to about detector.threshold ≈ 0.6 - 0.8.


Examples
========

#. Running magsdetector offline with a multiplexed miniseed volume and an
   inventory xml file. Neither a messaging nor a database connection is required.

   .. code-block:: sh

      magsdetector --inventory-db inventory.xml -I test-sorted.mseed --offline
