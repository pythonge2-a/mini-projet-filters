class Filters:
    def __init__(self):
        self.snk = None  # Définir snk comme `None` par défaut

    def get_snk(self):
        if self.snk is None:
            from .snk import SallenAndKey  # Import différé
            self.snk = SallenAndKey()
        return self.snk
