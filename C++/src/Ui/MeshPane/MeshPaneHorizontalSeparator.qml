import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

RowLayout {
    Layout.alignment: Qt.AlignHCenter
    Layout.fillWidth: true
    Layout.preferredHeight: 1.0
    Item {
        Layout.preferredHeight: 1
        Layout.fillWidth: true
    }
    Rectangle{
        Layout.preferredHeight: 2.0
        Layout.fillWidth: true
        radius: 20
        color: "gray"
    }
    Item {
        Layout.preferredHeight: 1
        Layout.fillWidth: true
    }
}
