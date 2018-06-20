class Base(object):

    def clone(self):
        """Copy object.

        Args:
            self (class): class name.

        Returns:
            copy (class): class inheriting self's attributes.

        """
        copy = self.__class__()
        for col in self.__table__.columns:
            val = getattr(self, col.name)
            setattr(copy, col.name, val)
        return copy

    def update_from_dict(self, d):
        """Add properties to a class.

        Args:
            self (class): class name.
            d (dict): attributes in a dictionary.

        Dynamically adds property from dictionary to class-attribute.
        
        """
        for k, v in d.items():
            setattr(self, k, v)

from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base(cls=Base)
