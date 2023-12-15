import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

SwipeView {
    id: view
    Layout.fillWidth: true
    Layout.preferredHeight: 200
    currentIndex: 0
    padding: 0

    TetWildControl{}
    MeshRenderingControl{}

    Rectangle {
        Layout.fillWidth: true
        Layout.fillHeight: true
        color: "red"
    }
}
