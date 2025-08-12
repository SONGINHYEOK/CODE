#!/bin/bash

echo "Setting up OMNIDRUG-AI Platform..."

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r backend/requirements.txt

# Copy environment file
cp .env.example .env
echo "Please edit .env file with your configuration"

# Create PostgreSQL database
createdb omnidrug
createuser omnidrug_user

# Run migrations
cd backend
python manage.py makemigrations core projects molecules screening
python manage.py migrate

# Create superuser
echo "Creating superuser..."
python manage.py createsuperuser

# Collect static files
python manage.py collectstatic --noinput

echo "Setup complete!"
echo "To start the development server:"
echo "  cd backend"
echo "  python manage.py runserver"
echo ""
echo "To start Celery worker:"
echo "  celery -A omnidrug worker -l info"
echo ""
echo "To start Flower (Celery monitoring):"
echo "  celery -A omnidrug flower"