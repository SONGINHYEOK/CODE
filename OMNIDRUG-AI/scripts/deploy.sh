#!/bin/bash

echo "Deploying OMNIDRUG-AI Platform..."

# Build Docker images
docker-compose -f docker/docker-compose.yml build

# Run database migrations
docker-compose -f docker/docker-compose.yml run backend python manage.py migrate

# Collect static files
docker-compose -f docker/docker-compose.yml run backend python manage.py collectstatic --noinput

# Start services
docker-compose -f docker/docker-compose.yml up -d

echo "Deployment complete!"
echo "Access the application at http://localhost"
echo "Access Flower (Celery monitoring) at http://localhost:5555"