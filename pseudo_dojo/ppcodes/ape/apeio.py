from __future__ import division, print_function

import os.path
import collections
import numpy as np

from pseudo_dojo.core import RadialFunction, RadialWaveFunction, plot_aepp, plot_logders

__version__ = "0.1"
__status__ = "Development"
__date__ = "$April 26, 2013M$"

# FIXME circular dep
class ApeError(Exception):
    """Base class for APE Exceptions"""


class ApeEvent(object):
    """
    Example (JSON)
    <ApeEvent> 
        class: "ScfConvergenceWarning",
        tolname: tolwfr,
        actual_tol: 1.0e-8,
        required_tol: 1.0e-10, 
        nstep: 50,
        message: "The human-readable message goes here!"
    </ApeEvent>
    """
    @staticmethod
    def from_string(string, lineno):
        """Constructs an APE event given a string and the line number."""
        d = json.loads(string)
        cls = d.pop["class"]
        assert "lineno" not in d
        d["lineno"] = lineno
        return cls(d)

    def __init__(self, **kwargs):
        self._kwargs = kwargs.copy()
        self.message = kwargs.pop("message")
        self.lineno  = kwargs.pop("lineno")
        self.data = kwargs

    def __str__(self):
        return "%s:\n%s" % (self.lineno, self.message)

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def baseclass(self):
        for cls in _BASE_CLASSES:
            if isinstance(self, cls):
                return cls
        raise ValueError("Cannot determine the base class of %s" % self.__class__.__name__)

    def iscritical(self):
        """
        True if event is critical namely that if this event should be analyzed in 
        more detail to understand what action should be performed
        """
        return False

##########################################################################################


class Comment(ApeEvent):
    """Base class for Comment events."""


class Error(ApeEvent):
    """Base class for Error events."""

    @property
    def iscritical(self):
        return True


class Bug(ApeEvent):
    """Base class for Bug events."""

    @property
    def iscritical(self):
        return True


class Warning(ApeEvent):
    """
    Base class for Warning events (the most important class).
    Developers should subclass inherit from this class to define the different exceptions 
    raised by the code and the possible actions that can be performed.
    """

    # FIXME: for the moment we tag a warning as critical, then, once we migrate to the
    # JSON-like format, only CriticalWarnings will trigger some kind of action.
    @property
    def iscritical(self):
        return True

# Register the concrete base classes.
_BASE_CLASSES = [
    Comment,
    Error,
    Bug,
    Warning,
]

##########################################################################################


class EventList(collections.Iterable):
    """Iterable storing the events produced by APE."""
    def __init__(self, filename, events=None):
        """
        Args:
            filename:
                Name of the file
            events:
                List of Event objects
        """
        self.filename = os.path.abspath(filename)
        self._events = []
        self._events_by_baseclass = collections.defaultdict(list)

        if events is not None:
            for e in events:
                self.append(event)

    def __len__(self):
        return len(self._events)

    def __iter__(self):
        return self._events.__iter__()

    def __str__(self):
        lines = [self.filename+":",]
        for event in self:
            lines.append(str(event))
        return "\n".join(lines)

    def append(self, event):
        """Add an event to the list."""
        self._events.append(event)
        self._events_by_baseclass[event.baseclass].append(event)

    @property
    def critical_events(self):
        """List of critical events."""
        return [e for e in self if e.iscritical]

    @property
    def comments(self):
        return self.select(Comment)

    @property
    def errors(self):
        return self.select(Error)

    @property
    def bugs(self):
        return self.select(Bug)

    @property
    def warnings(self):
        return self.select(Warning)

    def select(self, base_class, only_critical=False):
        """
        Return list of events that inherits from class base_class

        only_critical:
            if True, only critical events are returned.
        """
        if only_critical:
            return [e for e in self._events_by_baseclass[base_class] if e.iscritical]
        else:
            return self._events_by_baseclass[base_class][:]

##########################################################################################


class ApeOutputParser(object):
    """
    Parses the output file produced by APE, extracts the most
    important results and the list of significant events.
    """
    Error = ApeError

    def __init__(self, out_fname):
        """
        Args:
            out_fname:
                Path to the APE output file.
        """
        self.out_fname = os.path.abspath(out_fname)
        self.dirpath = os.path.abspath(os.path.dirname(self.out_fname))

        with open(self.out_fname, "r") as fh:
            self.out_lines = tuple(fh.readlines())

    @property
    def events(self):
        """List of APE events."""
        try:
            return self._events
        except AttributeError:
            self._events = ape_parse_events(self.out_fname)
            return self._events

##########################################################################################


def ape_parse_events(out_lines):
    """
    Read and parse the main output file of APE.

    Args:
        out_lines:
            List of string containing the output file.

    Returns:
        EventList instance.
    """
    START_TAG = "<ApeEvent>"
    STOP_TAG = "<\ApeEvent>"

    events = EventList(out_fname)

    in_event, s = False, None

    for (l, line) in enumerate(out_lines):
        if not in_event:
            if line.startswith(START_TAG):
                in_event = True
                if s is None:
                    # First event found.
                    s = ""
                else:
                    # Parse the previous string to generate the appropriate event.
                    events.append(ApeEvent.from_string(s, lineno))
                lineno = l
        else:
            if line.startswith(STOP_TAG):
                in_event = False
            else:
                # Add new line.
                s += line

    if s:
        events.append(ApeEvent.from_string(s, lineno))

    return events

##########################################################################################

#class ApePlotter(object):
#
#    def build_figs(self, dirpath, **kwargs):
#        """
#        Args:
#
#        ==============  ==============================================================
#        kwargs          Meaning
#        ==============  ==============================================================
#        show            True to show the figure (Default).
#        ==============  ==============================================================
#        """
#        if "show" not in kwargs:
#            kwargs["show"] = False
#
#        self.figs = figs = {}
#
#        figs["wave"] = ape_plot_waves(dirpath, **kwargs) 
#
#        figs["logd"] = ape_plot_logders(dirpath, **kwargs)
#
#        figs["pot"] = ape_plot_potentials(dirpath, **kwargs)
#
#        figs["den"] = ape_plot_densities(dirpath, **kwargs)

##########################################################################################


def ape_read_waves(dirpath):
    """Read the APE radial wavefunctions located in directory dirpath."""
    waves = {}
    for filename in os.listdir(dirpath):
        if not filename.startswith("wf-"): 
            continue

        path = os.path.join(dirpath, filename)
        tokens = filename.split("-")
        assert len(tokens) == 2
        state = tokens[1]
        # TODO check spin and spinor
        # Load [r, Radial wavefunctions and first derivative] in data.
        data = np.loadtxt(path)
        waves[state] = RadialWaveFunction(state, state, data[:,0], data[:,1])

    return waves


def ape_plot_waves(dirpath, **kwargs): 
    """
    Uses Matplotlib to plot the radial wavefunction (AE vs PP)

    Args:
        dirpath:

    ==============  ==============================================================
    kwargs          Meaning
    ==============  ==============================================================
    title           Title of the plot (Default: None).
    show            True to show the figure (Default: True).
    savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
    ==============  ==============================================================

    Returns:
        `matplotlib` figure.
    """
    dirs = os.listdir(os.path.abspath(dirpath))

    ae_waves, pp_waves = None, None

    if "ae" in dirs:
        ae_waves = ape_read_waves(os.path.join(dirpath, "ae"))

    if "pp" in dirs:
        pp_waves = ape_read_waves(os.path.join(dirpath, "pp"))

    if ae_waves is not None:
        fig = plot_aepp(ae_waves, pp_funcs=pp_waves, **kwargs)
    else:
        print("Cannot find AE waves in dirpath %s" % dirpath)
        fig = None

    return fig



def ape_read_potentials(dirpath):
    """Read the APE radial potentials located in directory dirpath."""
    pots = {}
    for filename in os.listdir(dirpath):
        if not filename.startswith("v_"):
            continue
        path = os.path.join(dirpath, filename)
                                                                          
        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
                                                                       
        pots[filename] = RadialFunction(filename, data[:,0], data[:,1])
    return pots


def ape_plot_potentials(dirpath, **kwargs): 
    """
    Uses Matplotlib to plot the potentials (AE vs PP)

    Args:
        dirpath:

    ==============  ==============================================================
    kwargs          Meaning
    ==============  ==============================================================
    title           Title of the plot (Default: None).
    show            True to show the figure (Default: True).
    savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
    ==============  ==============================================================

    Returns:
        `matplotlib` figure.
    """
    #raise NotImplementedError("")
    dirs = os.listdir(os.path.abspath(dirpath))

    ae_pots, pp_pots = None, None

    if "ae" in dirs:
        ae_pots = ape_read_potentials(os.path.join(dirpath, "ae"))

    if "pp" in dirs:
        pp_pots = ape_read_potentials(os.path.join(dirpath, "pp"))

    if ae_pots is not None:
        fig = plot_aepp(ae_pots, pp_funcs=pp_pots, **kwargs)
    else:
        print("Cannot find potentials in dirpath %s" % dirpath)
        fig = None

    return fig



def ape_read_densities(dirpath):
    """Read APE AE densities and tau located in directory dirpath."""
    dens = {}
    for filename in os.listdir(dirpath):
        if filename not in ["density", "tau",]:
            continue
        path = os.path.join(dirpath, filename)

        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
        dens[filename] = RadialFunction(filename, data[:,0], data[:,1])
        #print(filename, dens[filename].values)
    return dens


def ape_plot_densities(dirpath, **kwargs): 
    """
    Uses Matplotlib to plot the densities (AE vs PP)

    Args:
        dirpath:
            Directory path.

    ==============  ==============================================================
    kwargs          Meaning
    ==============  ==============================================================
    title           Title of the plot (Default: None).
    show            True to show the figure (Default: True).
    savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
    ==============  ==============================================================

    Returns:
        `matplotlib` figure.
    """
    dirs = os.listdir(os.path.abspath(dirpath))

    ae_dens, pp_dens = None, None

    if "ae" in dirs:
        ae_dens = ape_read_densities(os.path.join(dirpath, "ae"))

    if "pp" in dirs:
        pp_dens = ape_read_densities(os.path.join(dirpath, "pp"))

    if ae_dens is not None:
        fig = plot_aepp(ae_dens, pp_funcs=pp_dens, **kwargs)
    else:
        print("Cannot find densities in dirpath %s" % dirpath)
        fig = None

    return fig



def ape_read_logders(dirpath):
    """Reads the APE log derivatives located in directory dirpath."""
    ae_logders, pp_logders = {}, {}
    for filename in os.listdir(dirpath):
        if not filename.startswith("ld-"):
            continue
        path = os.path.join(dirpath, filename)

        tokens = filename.split("-")
        assert len(tokens) == 2
        l_name = tokens[1]

        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
        ae_logders[l_name] = RadialFunction(l_name, data[:,0], data[:,1])
        pp_logders[l_name] = RadialFunction(l_name, data[:,0], data[:,2])

    return ae_logders, pp_logders


def ape_plot_logders(dirpath, **kwargs): 
    """
    Uses Matplotlib to plot the radial wavefunction (AE vs PP)

    Args:
        dirpath:
            Directory path.

    ==============  ==============================================================
    kwargs          Meaning
    ==============  ==============================================================
    title           Title of the plot (Default: None).
    show            True to show the figure (Default: True).
    savefig:        'abc.png' or 'abc.eps'* to save the figure to a file.
    ==============  ==============================================================

    Returns:
        `matplotlib` figure.
    """
    dirpath = os.path.abspath(dirpath)
    if "tests" in os.listdir(dirpath):
        dirpath = os.path.join(dirpath, "tests")

    ae_logders, pp_logders = ape_read_logders(dirpath)

    fig = plot_logders(ae_logders, pp_logders, **kwargs)
    return fig



def ape_check_ghosts(out_lines):
    """
    Check for the presence of ghost states. Reads data from ape.out.

    Example::

          Ghost state analysis:
            State: 3s
              KB energy < 0; Eref < E0       =>  No ghost states
              Local potential eigenvalues:   -0.1182 (E0)     0.0000 (E1)
              Reference energy:              -0.3981 (Eref)
            State: 3d
              KB energy < 0; Eref = E0 = 0   =>  Unable to determine
              Local potential eigenvalues:    0.0000 (E0)     0.0000 (E1)
              Reference energy:               0.0000 (Eref)
            State: 4f
              KB energy < 0; Eref > E0       =>  Ghost state found
              Local potential eigenvalues:   -3.0076 (E0)     0.0000 (E1)
              Reference energy:               0.0000 (Eref)

            Localization radii [b]:
    """
    out_lines = list(out_lines[:])

    SENTINEL = "Ghost state analysis:"
    for (i, line) in enumerate(out_lines):
        if line.strip() == SENTINEL:
            out_lines = out_lines[i+1:]
            break
    else:
        raise ApeError("Cannot find sentinel %s in APE output" % SENTINEL)

    ghosts = {}
    while True:
        line = out_lines.pop(0).strip()
        if not line: 
            break

        if line.startswith("State:"):
            state = line.split("State:")[1].strip()
            line = out_lines.pop(0).strip()
            result = line.split("=>")[-1].lower().strip()
            result = result.split()[0]
            # See states.f90
            assert result in ["no", "ghost", "unable", "illdefined"]
            if result == "ghost":
                ghosts[state] = result

    return ghosts


def ape_check_ppeigen(out_lines):
    """
    Check the quality of the pseudo eigenvalues.

    Pseudopotentials Self-Consistency:
      State  Eigenvalue [H ]    Norm Test   Slope Test
        3s        -0.39812      1.0000003   0.9999983
        3p        -0.15331      0.9999993   0.9999941
        3d         0.00000      1.0033912   1.0017468
        4f         0.00000      0.9912183   0.9956052

    """
    class PseudoScfTest(collections.namedtuple("PseudoScfTest", "state, eig, norm, slope")):
        def __new__(cls, **kwargs):
            # Type conversion
            for k, v in kwargs.items():
                if k == "state":
                    kwargs[k] = str(v)
                else:
                    kwargs[k] = float(v)
            return super(PseudoScfTest, cls).__new__(cls, **kwargs)

    out_lines = list(out_lines[:])

    SENTINEL = "Pseudopotentials Self-Consistency:"
    for (i, line) in enumerate(out_lines):
        if line.strip() == SENTINEL:
            out_lines = out_lines[i+2:]
            break
    else:
        raise ApeError("Cannot find sentinel %s in APE output" % SENTINEL)

    scf_tests = {}
    while True:
        line = out_lines.pop(0).strip()
        if not line: 
            break

        state, eig, norm, slope = line.split()

        scf_tests[state] = PseudoScfTest(state=state, eig=eig, norm=norm, slope=slope)


def ape_read_dipoles(out_lines):
    """
    Args:
        out_lines:
            List of strings with the APE output.
    Returns:
        List of Dipole objects.
    """
    # Dipole Matrix Elements:
    #      States           AE         PS
    #     3s -- 3p        2.2998     2.3074
    #     3s -- 3d        0.0003     0.0003
    #     3s -- 4f        0.0000     0.0000
    #     3p -- 3d        0.0008     0.0008
    #     3p -- 4f        0.0000     0.0000
    #     3d -- 4f       98.7760    98.7760
    out_lines = list(out_lines[:])
                                                                           
    SENTINEL = "Dipole Matrix Elements:"
    for (i, line) in enumerate(out_lines):
        if line.strip() == SENTINEL:
            out_lines = out_lines[i+2:]
            break
    else:
        raise ApeError("Cannot find sentinel %s in APE output" % SENTINEL)

    dipoles = {}
    while True:
        line = out_lines.pop(0).strip()
        if not line: 
            break
        tokens = line.split()

        ae_r, ps_r = map(float, [tokens[-2], tokens[-1]])
        trans = " ".join(tokens[:3])
        dipoles[trans] = {"AE": ae_r, "PS": ps_r, "AE-PS": (ae_r-ps_r)}

        #Dipole(istate, ostate, aeres, ppres):

    return dipoles

