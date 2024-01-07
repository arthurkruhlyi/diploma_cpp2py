class Material:
    def __init__(self, melting_temperature, boiling_temperature, eps):
        self.melting_temperature = melting_temperature
        self.boiling_temperature = boiling_temperature
        self.eps = eps

"""

{
public:
	Matherial(double, double, double);
	~Matherial();

	void SetTemperatureOfMelting(double);
	double GetTemperatureOfMelting ();

	void SetBoilingTemperature(double);
	double GetBoilingTemperature();

	void SetEps(double);
	double GetEps();

private:
	double TemperatureOfMelting;
	double BoilingTemperature;
	double Eps;
};

class Observable:
    def __init__(self):
        self.observers = []

    def register(self, observer):
        self.observers.append(observer)

    def unregister(self, observer):
        self.observers.remove(observer)

    def notify(self, f):
        def wrapper(*args, **kwargs):
            f(*args, **kwargs)
            for observer in self.observers:
                observer()

        return wrapper


o = Observable()
"""
