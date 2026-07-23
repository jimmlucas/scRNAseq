#!/bin/bash
set -e

# Use provided env vars or fallback defaults
APP_USER=${USER:-rstudio}
APP_USERID=${USERID:-1000}
APP_GROUPID=${GROUPID:-1000}
PASSWORD=${PASSWORD:-rstudio}
HOME_DIR="/home/${APP_USER}"

# Create group if it doesn't exist
if ! getent group ${APP_GROUPID} >/dev/null; then
    groupadd -g ${APP_GROUPID} ${APP_USER}
fi

# Create user if it doesn't exist
if ! id -u ${APP_USER} >/dev/null 2>&1; then
    useradd -m -s /bin/bash -u ${APP_USERID} -g ${APP_GROUPID} ${APP_USER}
fi

# Set user password
echo "${APP_USER}:${PASSWORD}" | chpasswd

# Set permissions
mkdir -p ${HOME_DIR}
chown -R ${APP_USER}:${APP_GROUPID} ${HOME_DIR}

# Start RStudio server as that user
exec /init